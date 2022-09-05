non_na_mean(vv) = mean(vv[.!(isnan.(vv))])

non_na_mean(vv) = mean(vv[.!(isnan.(vv))])

function fix_metadata_colnames!(longdata_df::DataFrame)

    ## This function was created to fix small inconsistencies and typos on the longdata table.
    ## As original data cleanup and wrangling progresses, it will be altered and ultimately deprecated.
    ## Comments are provided to help review transformations

    # 1. Rename "[Collinsella]_massiliensis" to "Collinsella_massiliensis"
    longdata_df[longdata_df.variable .== "[Collinsella]_massiliensis", :variable] .= "Collinsella_massiliensis"

    return(longdata_df)
    
end


function check_longdata_metaduplicates!(longdata_df::DataFrame; remove_duplicates=true)
    n_unique_rows = size(unique(longdata_df[:, 1:3]),1)    
    
    if (n_unique_rows != size(longdata_df, 1))
        n_nonunique = size(longdata_df, 1) - n_unique_rows

        if remove_duplicates    
            @warn "Long Dataframe contains non-unique rows! Argument `remove_duplicates` set to `true`. Removing duplicated data."
            @warn "After removal, $(n_unique_rows) will remain. $(n_nonunique) rows were rmoved from the original $(size(longdata_df, 1))"
            rawDf = unique!(longdata_df)
            return(longdata_df)
        else
            @warn "Long Dataframe contains non-unique rows! Argument `remove_duplicates` set to `false`. Returning source data."
            return(longdata_df)
        end #end if remove_duplicates
    end # end if n_unique_rows
end # end function

function build_metadata_prediction_df(base_df, inputs::Vector{Symbol}, targets::Vector{Symbol})

    subjects = unique(base_df.subject)
    has_stool = findall(.!(ismissing.(base_df.sample)))
    has_cogscore = findall(.!(ismissing.(base_df.cogScore)))

    subjects_stool_idx = [ intersect(findall(base_df.subject .== el), has_stool) for el in subjects ]
    subjects_cogscore_idx = [ intersect(findall(base_df.subject .== el), has_cogscore) for el in subjects ]

    line_combinations = vcat( [ vec(collect(Base.product(subjects_stool_idx[i], subjects_cogscore_idx[i]))) for i in 1:length(subjects) ]... )
    line_combinations = line_combinations[ [ el[2] > el[1] for el in line_combinations ] ]

    targets_df = @chain base_df[ [el[2] for el in line_combinations], :] begin
        select( [:subject, :timepoint, :ageMonths, targets...] )
        rename!( [:timepoint => :futureTimepoint, :ageMonths => :futureAgeMonths, (targets .=> Symbol.("future" .* uppercasefirst.(String.(targets))))...] )
    end

    inputs_df = @chain base_df[ [el[1] for el in line_combinations], :] begin
        select( Not(:subject) )
    end

    prediction_df = hcat(targets_df, inputs_df)
    prediction_df.ageMonthsDelta = prediction_df.futureAgeMonths .- prediction_df.ageMonths 
    prediction_df.timepointDelta = prediction_df.futureTimepoint .- prediction_df.timepoint 
    
    select!(prediction_df, [:subject, :ageMonths, :timepoint, :futureAgeMonths, :futureTimepoint, :ageMonthsDelta, :timepointDelta, targets..., Symbol.("future" .* uppercasefirst.(String.(targets)))..., inputs...])
    sort!(prediction_df, [:subject, :timepoint, :futureTimepoint])
    
    return prediction_df

end

function tryparsecol(T, col)
    newcol = tryparse.(T, col)
    return any(isnothing, newcol) ? col : newcol
end

function compute_tietjenmoore(data, k, n)
    r_all = abs.(data .- mean(data))
    filteredData = data[sortperm(r_all)][1:(n-k)]
    ksub = (filteredData .- mean(filteredData))

    return( sum(ksub .^ 2) / sum(r_all .^ 2) )
end

function test_tietjenmoore(dataSeries,k, n, alpha, sim_trials)
    ek = compute_tietjenmoore(dataSeries,k, n)
    simulation = [ compute_tietjenmoore(randn(length(dataSeries)), k, n) for i in 1:sim_trials ]
    Talpha=quantile(simulation,alpha)

    return(ek, Talpha)
end

function univariate_tietjenmoore(values::Vector{T} where T <: Real, k::Int64; alpha = 0.05, sim_trials = 100000)

    @info "----- Begin Tietjen-Moore Outlier test -----"
    @info "H0: There are no outliers in the data set."
    @info "H1: There are exactly k outliers in the data set\n"

    n = length(values)

    L_set, L_critical = test_tietjenmoore(values, k, n, alpha, sim_trials)
    @info "Set L-statistic for $n samples and $k outliers with mode $(mode): $(round(L_set, digits = 4))"
    @info "Critical L for $n samples and $k outliers with mode $(mode): $(round(L_critical, digits = 4))\n"

    if L_set < L_critical
        @info "L_set < L_critical !"
        @info "**SUCCESSFUL REJECTION OF H0** with confidence level $alpha" 
        r_all = abs.(values .- mean(values))
        outlier_indexes = sortperm(r_all)[(n-k+1):end]
        return outlier_indexes
    else
        @warn "L_set > L_critical !"
        @warn "**CANNOT REJECT H0** with confidence level $alpha"
        return Int64[]
    end # endif L_set < L_critical
end # end function

function try_outliers(f, data, n)
    for i in n:-1:1
        outlier_idx = f(data, i)
        if length(outlier_idx) != 0
            return(i, outlier_idx)
        end
    end
    return 0, Int64[]
end

function unstack_techreplicates(raw_data, retain_first=true)

    unstackedDf = retain_first ?
    begin 
        @chain raw_data begin
            groupby([:subject, :timepoint, :variable])
            combine(:value=>first)
            unstack(:variable, :value_first)
        end #end @chain
    end : #end if_true
    unstack(raw_data, :variable, :value; allowduplicates=true)

    return unstackedDf
end #end function

function report_classification_merit(classification_results::Dict)
    classification_merits = DataFrame(
    :Trial => collect(1:classification_results[:n_trials]),
    :Train_ACC => classification_results[:train_accuracies],
    :Test_ACC => classification_results[:test_accuracies],
    )

    push!(classification_merits, Dict(
        :Trial => 0,
        :Train_ACC => mean(classification_merits.Train_ACC),
        :Test_ACC => mean(classification_merits.Test_ACC)
    ))

    return classification_merits
end

function report_regression_merit(regression_results::Dict)
    regression_merits = DataFrame(
        :Trial => collect(1:regression_results[:n_trials]),
        :Train_MAE => regression_results[:train_maes],
        :Test_MAE => regression_results[:test_maes],
        :Train_MAPE => regression_results[:train_mapes],
        :Test_MAPE => regression_results[:test_mapes],
        :Train_COR => regression_results[:train_correlations],
        :Test_COR => regression_results[:test_correlations]
        )

    push!(regression_merits, Dict(
        :Trial => 0,
        :Train_MAE => mean(regression_merits.Train_MAE),
        :Test_MAE => mean(regression_merits.Test_MAE),
        :Train_MAPE => mean(regression_merits.Train_MAPE),
        :Test_MAPE => mean(regression_merits.Test_MAPE),
        :Train_COR => mean(regression_merits.Train_COR),
        :Test_COR => mean(regression_merits.Test_COR)
    ))
    
    return regression_merits
end

function build_confusion_matrix(classification_results::Dict, trial::Int)

    y = classification_results[:inputs_outputs][2][classification_results[:dataset_partitions][trial][2]]
    yhat = MLJ.predict_mode(classification_results[:models][trial], classification_results[:inputs_outputs][1][classification_results[:dataset_partitions][trial][2],:])
    confmat = ConfusionMatrix()(yhat, y)

    return confmat
end

function average_confusion_matrix(classification_results::Dict)

    matrix_vector = [ build_confusion_matrix(classification_results, i).mat for i in 1:classification_results[:n_trials] ]
    concatenated_matrix = vcat([ vec(mat)' for mat in matrix_vector ]...) # Column order is TN, FP, FN, TP !
    averages = [ sum(concatenated_matrix[:,i])/sum(concatenated_matrix) for i in 1:4]

    return averages
end

function confmatrix2barplot(classification_results::Dict)

    confmat_inputs = (
        x = [1, 2, 2, 1], # Column order is TN, FP, FN, TP !
        value = average_confusion_matrix(classification_results),
        grp = [1, 2, 2, 1],
        color = [1, 2, 3, 4]
    )

    return confmat_inputs
end

function regression_bestprediction(regression_results::Dict)

    selected_trial = regression_results[:selected_trial][2]
    X, y = regression_results[:inputs_outputs]
    train, test = regression_results[:dataset_partitions][selected_trial]
    slope_correction = regression_results[:slope_corrections][selected_trial]
    selected_machine = regression_results[:models][selected_trial]
    yhat = GLM.predict(slope_correction, DataFrame( :yhat => MLJ.predict(selected_machine, X)))

    return y, yhat, train, test
end