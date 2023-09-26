function runlms(comm, outfile;
                modcols=[ "cogScore", "ageMonths", "read_depth", "education"],
                formula=@formula(bug ~ cogScore + ageMonths + read_depth + education),
                keyfeature=1,
                prevalence_threshold = 0.2,
                model_kind = :regular # or :logistic
    )
    indf = DataFrame(get(comm))[:, modcols]

    lmresults = DataFrame(ThreadsX.map(features(comm)) do feature
        @debug feature
        default_nt = (; Name = modcols[keyfeature], feature = name(feature), coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN)            

        ab = vec(abundances(comm[feature, :]))
        over0 = ab .> 0
        sum(over0) / size(indf, 1) > prevalence_threshold || return default_nt
 
        # ab = collect(indf[!, feature] .+ (minimum(indf[over0, feature])) / 2) # add half-minimum non-zerovalue

        df = indf[!, modcols]
        if model_kind == :regular
            df.bug = asin.(sqrt.(ab))
            lmod = lm(formula, df; dropcollinear=false)
        elseif model_kind == :logistic
            df.bug = over0
            lmod = glm(formula, df, Binomial(), ProbitLink(); dropcollinear=false)
        else
            throw(ArgumentError("No model_kind: $model_kind"))    
        end
        ct = DataFrame(coeftable(lmod))
        ct.feature .= name(feature)
        rename!(ct, "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
        model_kind == :regular && rename!(ct, "Pr(>|t|)"=>"pvalue", "t"=>"stat")
        model_kind == :logistic && rename!(ct, "Pr(>|z|)"=>"pvalue", "z"=>"stat")
        select!(ct, Cols(:feature, :Name, :))

        return NamedTuple(only(filter(row-> row.Name == modcols[keyfeature], eachrow(ct))))
    end)

    subset!(lmresults, "pvalue"=> ByRow(!isnan))
    DataFrames.transform!(lmresults, :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(lmresults, :qvalue)

    CSV.write(outfile, lmresults)
    lmresults
end

