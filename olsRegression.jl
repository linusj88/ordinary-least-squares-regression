using Base, Distributions, DataFrames, RDatasets, Plots, CSV, LinearAlgebra

lm = function(data::DataFrame, dependent, independent, significans_level = 0.05)

    # Create design matrix
    d = data
    n = size(d)[1]

    # Create design matrix
    X = convert(Array, [ones(n) d[independent]])
    y = d[dependent]
    parameter = vcat(:intercept, independent)

    # Fit model by solving normal equations.
    Xt = X'
    invXtX = inv(Xt * X)
    H = X * invXtX * Xt
    pred = H * y
    param = invXtX * Xt * y
    p = length(param)

    # Compute statistics of the fit
    res = (y .- pred) .^2
    rss = sum(res)
    df = n - p
    mse = rss / df
    seParam = (diag(invXtX .* mse)) .^.5

    # Perform t-test for the significans of the parameters
    z = param ./ seParam
    t = TDist(df)
    pvalue = ones(p)
    for i in 1:p
        if z[i] < 0
            pvalue[i] = cdf(t, z[i])
        else
            pvalue[i] = 1 - cdf(t, z[i])
        end
    end

    fit = DataFrame(Parameter = parameter, Value = param, Std_error = seParam, pValue = pvalue)

    # Compute statistics about the fit
    sstotal = sum((y .- mean(y)) .^2)
    rsquared = 1 - rss / sstotal
    rsquaredadj = 1 - (1 - rsquared) * (n - 1) / df
    variation = DataFrame(Meassure = ["R²", "R² adjusted", "√MSE", "RSS", "Degrees of freedom"], Value = [rsquared, rsquaredadj, sqrt(mse), rss, df])

    return (fit = fit, variation = variation, residuals = res, predictions = pred)
end


# Load data
d = CSV.read("C:/Users/oskar/Documents/challenger.csv", delim = "\t")
d.freq = d.rings_under_stress ./ 6
print(sort!(d, cols = [:fahrenheit]))

# Define dependent and independent variables and variables to not include
dependent = :freq
independent = [:fahrenheit]
println()
m = lm(d, dependent, independent)

#
plot(d.fahrenheit,
        d.freq,
        seriestype = :scatter,
        xlabel = "°F",
        xlims = (20,83),
        xticks = 20:5:82,
        ylims = (-.05,1.05),
        yticks = 0:0.1:1,
        ylabel = "O-rings",
        markersize = 3,
        markercolor = :black)
n = (m.fit).Value
x = 20:.01:81
y = collect(n[1] .+ n[2] .* x)
plot!(x, y, seriestype = :line, color = :red)

# Fit a polynomial
d.fahrenheit2 = d.fahrenheit .^2
m2 = lm(d, :freq, [:fahrenheit, :fahrenheit2])
y2 = m2.fit.Value[1] .+ m2.fit.Value[2] .* x .+ m2.fit.Value[3] .* x .^2
plot!(x, y2)

scatter(m2.predictions, m2.residuals)

i = d.fahrenheit .<= 70
m3 = lm(d[.! i,:], :freq, :fahrenheit)
scatter(m3.predictions, m3.residuals)
