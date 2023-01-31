function runlms(indf, outfile, featurecols;
                modelcols=[ "cogScore", "ageMonths", "read_depth", "education"],
                formula=@formula(bug ~ cogScore + ageMonths + read_depth + education),
                keyfeature=1)

    lmresults = DataFrame(ThreadsX.map(featurecols) do feature
        @debug feature
            
        over0 = indf[!, feature] .> 0
        sum(over0) / size(indf, 1) > 0.20 || return (Name = modelcols[keyfeature], feature, coef = NaN, std_err = NaN, t = NaN, pvalue = NaN, lower_95=NaN, upper_95=NaN)
        # ab = collect(indf[!, feature] .+ (minimum(indf[over0, feature])) / 2) # add half-minimum non-zerovalue

        df = indf[:, modelcols]
        df.bug = asin.(sqrt.(indf[!, feature]))

        mod = lm(formula, df; dropcollinear=false)
        ct = DataFrame(coeftable(mod))
        ct.feature .= feature
        rename!(ct, "Pr(>|t|)"=>"pvalue", "Lower 95%"=> "lower_95", "Upper 95%"=> "upper_95", "Coef."=> "coef", "Std. Error"=>"std_err")
        select!(ct, Cols(:feature, :Name, :))
        return NamedTuple(only(filter(row-> row.Name == modelcols[keyfeature], eachrow(ct))))
    end)

    subset!(lmresults, "pvalue"=> ByRow(!isnan))
    DataFrames.transform!(lmresults, :pvalue => (col-> MultipleTesting.adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    sort!(lmresults, :qvalue)

    CSV.write(outfile, lmresults)
    lmresults
end