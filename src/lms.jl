function runlms(indf, outfile, featurecols;
                modelcols=[ "cogScore", "ageMonths", "read_depth", "education"],
                formula=@formula(bug ~ cogScore + ageMonths + read_depth + education),
                keyfeature=1,
                prevalence_threshold = 0.2,
                model_kind = :regular # or :logistic
    )

    lmresults = DataFrame(ThreadsX.map(featurecols) do feature
        @debug feature
        default_nt = (; Name = modelcols[keyfeature], feature, coef = NaN, std_err = NaN, stat = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN)            

        over0 = indf[!, feature] .> 0
        sum(over0) / size(indf, 1) > prevalence_threshold || return default_nt
 
        # ab = collect(indf[!, feature] .+ (minimum(indf[over0, feature])) / 2) # add half-minimum non-zerovalue

        df = indf[:, modelcols]
        if model_kind == :regular
            df.bug = asin.(sqrt.(indf[!, feature]))
            lmod = lm(formula, df; dropcollinear=false)
        elseif model_kind = :logistic
            df.bug = over0
            lmod = glm(formula, df, Binomial(), ProbitLink(); dropcollinear=false)
        end
        ct = DataFrame(coeftable(lmod))
        ct.feature .= feature
        rename!(ct, "Pr(>|t|)"=>"pvalue", "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
        model_kind == :regular && rename!(ct, "t"=>"stat")
        model_kind == :logistic && rename!(ct, "z"=>"stat")
        select!(ct, Cols(:feature, :Name, :))

        return NamedTuple(only(filter(row-> row.Name == modelcols[keyfeature], eachrow(ct))))
    end)

    subset!(lmresults, "pvalue"=> ByRow(!isnan))
    DataFrames.transform!(lmresults, :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(lmresults, :qvalue)

    CSV.write(outfile, lmresults)
    lmresults
end

