using RDatasets, DataFrames, LinearAlgebra, Distributions, CSV

function interaction(df::DataFrame, a::Symbol, b::Symbol, interactionSymbol::Symbol)
        df2 = df
        df2[interactionSymbol] = df[a] .* df[b]
        return df2
end

function interaction!(df::DataFrame, a::Symbol, b::Symbol, interactionSymbol::Symbol)
        df[interactionSymbol] = df[a] .* df[b]
end

function categorize!(df, variable::Symbol)
        symbols = levels(df[variable])
        n = size(df)[1]
        for i in 1:length(symbols)-1
                j = df[variable] .== symbols[i]
                s = Symbol(symbols[i])
                df[s] = zeros(n::Int)
                df[j,s] = Int(1)
        end
end
