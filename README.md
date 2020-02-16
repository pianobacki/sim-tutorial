---
title: "Simulation Tutorial"
author: "Lisa DeBruine"
date: 2020-02-17
---

## Setup 

### Julia

Load the packages we'll be using in Julia. In pkg run the following to get the versions we're using:

* `add MixedModels#master`
* `add https://github.com/RePsychLing/MixedModelsSim.jl#master`

~~~~{.julia}
julia> using MixedModels        # run mixed models

julia> using MixedModelsSim     # simulation functions for mixed models

julia> using RCall              # call R functions from inside Julia

julia> using DataFrames, Tables # work with data tables

julia> using Random             # random number generator

julia> using CSV                # write CSV files

~~~~~~~~~~~~~





### R

Also load any packages we'll be using in R through `RCall()`.

~~~~{.julia}
R"""
library(ggplot2) # for visualisation
library(dplyr)   # for data wrangling
library(tidyr)   # for data wrangling
""";
~~~~~~~~~~~~~





### Define Custom functions

It's useful to be able to weave your file quickly while you're debugging, 
so set the number of simulations to a relatively low number while you're 
setting up your script and change it to a larger number when everything
is debugged.

~~~~{.julia}
nsims = 1000 # set to a low number for test, high for production
~~~~~~~~~~~~~


~~~~
1000
~~~~





#### Define: ggplot_betas

This function plots the beta values returned from `simulate_waldtests` using ggplot in R.
If you set a figname, it will save the plot to the specified file.

~~~~{.julia}
function ggplot_betas(sim, figname = 0, width = 7, height = 5) 

    beta_df = DataFrame(columntable(sim).β)

    R"""
        p <- $beta_df %>%
            gather(var, val, 1:ncol(.)) %>%
            ggplot(aes(val, color = var)) +
            geom_density(show.legend = FALSE) +
            facet_wrap(~var, scales = "free")

        if (is.character($figname)) {
            ggsave($figname, p, width = $width, height = $height)
        }

        p
    """
end
~~~~~~~~~~~~~


~~~~
ggplot_betas (generic function with 4 methods)
~~~~





## Existing Data

Load existing data from this morning's tutorial. Set the contrasts and run model 4 from the tutorial.

~~~~{.julia}
# load data
kb07 = MixedModels.dataset("kb07");

# set contrasts
contrasts = Dict(:spkr => HelmertCoding(), 
                 :prec => HelmertCoding(), 
                 :load => HelmertCoding());

# define formula
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );

# fit model
kb07_m = fit(MixedModel, kb07_f, kb07, contrasts=contrasts)
~~~~~~~~~~~~~


~~~~
Linear mixed model fit by maximum likelihood
 rt_trunc ~ 1 + spkr + prec + load + (1 | subj) + (1 + prec | item)
     logLik        -2 logLik          AIC             BIC       
 -1.43319251×10⁴  2.86638501×10⁴  2.86818501×10⁴  2.87312548×10⁴

Variance components:
             Column      Variance   Std.Dev.   Corr.
item     (Intercept)     133015.240 364.71254
         prec: maintain   63766.936 252.52116 -0.70
subj     (Intercept)      88819.437 298.02590
Residual                 462443.388 680.03190
 Number of obs: 1789; levels of grouping factors: 32, 56

  Fixed-effects parameters:
──────────────────────────────────────────────────────
                 Estimate  Std.Error  z value  P(>|z|)
──────────────────────────────────────────────────────
(Intercept)     2181.85      77.4681    28.16   <1e-99
spkr: old         67.879     16.0785     4.22   <1e-4 
prec: maintain  -333.791     47.4472    -7.03   <1e-11
load: yes         78.5904    16.0785     4.89   <1e-5 
──────────────────────────────────────────────────────
~~~~





### Simulate data with same parameters

Use the `simulate_waldtests()` function to run 1000 iterations of data sampled using the parameters from `m4`. Set up a random seed to make the simulation reproducible. You can use your favourite number.

To use multithreading, you need to set the number of cores you want to use. In Visual Studio Code, open the settings (gear icon in the lower left corner or cmd-,) and search for "thread". Set `julia.NumThreads` to the number of cores you want to use (at least 1 less than your total number).

~~~~{.julia}
# set seed for reproducibility
rng = MersenneTwister(8675309);

# run nsims iterations
kb07_sim = simulate_waldtests(rng, nsims, kb07_m, use_threads = true);
~~~~~~~~~~~~~





**Try**: Run the code above with and without `use_threads`.

Save all data to a csv file.

~~~~{.julia}
kb07_sim_df = sim_to_df(kb07_sim)

CSV.write("sim/kb07_sim.csv", kb07_sim_df)
~~~~~~~~~~~~~


~~~~
"sim/kb07_sim.csv"
~~~~





Plot betas in ggplot. In the code editor or Jupyter notebooks, you can omit the file name to just display the figure in an external window. 

~~~~{.julia}
# just display the image
# ggplot_betas(kb07_sim) 

# save the image to a file and display (display doesn't work in weave)
ggplot_betas(kb07_sim, "fig/kb07_betas.png");
~~~~~~~~~~~~~





In documents you want to weave, save the image to a file and use markdown to display the file. Add a semicolon to the end of the function to suppress creating the images in new windows during weaving.

![](fig/kb07_betas.png)


### Power calculation

The function `power_table()` from `MixedModelsSim` takes the output of `simulate_waldtests()` and calculates the proportion of simulations where the p-value is less than alpha for each coefficient. You can set the `alpha` argument to change the default value of 0.05 (justify your alpha ;).

~~~~{.julia}
power_table(kb07_sim)
~~~~~~~~~~~~~


~~~~
4×2 DataFrames.DataFrame
│ Row │ coefname       │ power   │
│     │ Symbol         │ Float64 │
├─────┼────────────────┼─────────┤
│ 1   │ (Intercept)    │ 1.0     │
│ 2   │ spkr: old      │ 0.991   │
│ 3   │ prec: maintain │ 1.0     │
│ 4   │ load: yes      │ 0.999   │
~~~~





### Change parameters

Let's say we want to check our power to detect effects of spkr, prec, and load 
that are half the size of our pilot data. We can set a new vector of beta values 
with the `β` argument to `simulate_waldtests`.

~~~~{.julia}
newβ = kb07_m.β
newβ[2:4] = kb07_m.β[2:4]/2

kb07_sim_half = simulate_waldtests(rng, nsims, kb07_m, β = newβ, use_threads = true);

power_table(kb07_sim_half)
~~~~~~~~~~~~~


~~~~
4×2 DataFrames.DataFrame
│ Row │ coefname       │ power   │
│     │ Symbol         │ Float64 │
├─────┼────────────────┼─────────┤
│ 1   │ (Intercept)    │ 1.0     │
│ 2   │ spkr: old      │ 0.529   │
│ 3   │ prec: maintain │ 0.928   │
│ 4   │ load: yes      │ 0.692   │
~~~~






# Simulating Data from Scratch


## simdat_crossed

The `simdat_crossed()` function from `MixedModelsSim` lets you set up a data frame with a specified experimental design. For now, it only makes fully balanced crossed designs, but you can generate an unbalanced design by simulating data for the largest cell and deleting extra rows. 

We will set a design where `subj_n` subjects per `age` group (O or Y) respond to `item_n` items in each of two `condition`s (A or B).

Your factors need to be specified separately for between-subject, between-item, and within-subject/item factors using `Dict` with the name of each factor as the keys and vectors with the names of the levels as values.

~~~~{.julia}
# put between-subject factors in a Dict
subj_btwn = Dict("age" => ["O", "Y"])

# there are no between-item factors in this design so you can omit it or set it to nothing
item_btwn = nothing

# put within-subject/item factors in a Dict
both_win = Dict("condition" => ["A", "B"])

# simulate data
dat = simdat_crossed(10, 30, 
                     subj_btwn = subj_btwn, 
                     item_btwn = item_btwn, 
                     both_win = both_win);
~~~~~~~~~~~~~






## Fit a model

Now you need to fit a model to your simulated data. Because the `dv` is just random numbers from N(0,1), there will be basically no subject or item random variance, residual variance will be near 1.0, and the estimates for all effects should be small. Don't worry, we'll specify fixed and random effects directly in `simulate_waldtests`. 

~~~~{.julia}
# set contrasts
contrasts = Dict(:age => HelmertCoding(), 
                 :condition => HelmertCoding());

f1 = @formula dv ~ 1 + age * condition + (1|item) + (1|subj);
m1 = fit(MixedModel, f1, dat, contrasts=contrasts)
~~~~~~~~~~~~~


~~~~
Linear mixed model fit by maximum likelihood
 dv ~ 1 + age + condition + age & condition + (1 | item) + (1 | subj)
   logLik   -2 logLik     AIC        BIC    
 -1732.2943  3464.5885  3478.5885  3514.2191

Variance components:
            Column     Variance   Std.Dev.  
item     (Intercept)  0.006272732 0.07920058
subj     (Intercept)  0.000000000 0.00000000
Residual              1.044877433 1.02219246
 Number of obs: 1200; levels of grouping factors: 30, 20

  Fixed-effects parameters:
───────────────────────────────────────────────────────────────
                          Estimate  Std.Error  z value  P(>|z|)
───────────────────────────────────────────────────────────────
(Intercept)            -0.0240362   0.0328606    -0.73   0.4645
age: Y                 -0.0142197   0.0295082    -0.48   0.6299
condition: B           -0.00221285  0.0295082    -0.07   0.9402
age: Y & condition: B  -0.00459405  0.0295082    -0.16   0.8763
───────────────────────────────────────────────────────────────
~~~~





## Simulate

Set a seed for reproducibility and specify β, σ, and θ.

~~~~{.julia}
rng = MersenneTwister(8675309);

new_beta = [0, 0.25, 0.25, 0]
new_sigma = 2.0
new_theta = [1.0, 1.0]

sim1 = simulate_waldtests(rng, nsims, m1, 
                        β = new_beta, 
                        σ = new_sigma, 
                        θ = new_theta,
                        use_threads = true);
~~~~~~~~~~~~~






## Explore simulation output


~~~~{.julia}
ggplot_betas(sim1, "fig/simbetas.png");
~~~~~~~~~~~~~





![](fig/simbetas.png)


## Power

~~~~{.julia}
power_table(sim1)
~~~~~~~~~~~~~


~~~~
4×2 DataFrames.DataFrame
│ Row │ coefname              │ power   │
│     │ Symbol                │ Float64 │
├─────┼───────────────────────┼─────────┤
│ 1   │ (Intercept)           │ 0.062   │
│ 2   │ age: Y                │ 0.132   │
│ 3   │ condition: B          │ 0.989   │
│ 4   │ age: Y & condition: B │ 0.056   │
~~~~






## Write a function to vary something

~~~~{.julia}
function mysim(subj_n, item_n, nsims = 1000, 
               beta  = [0, 0, 0, 0],
               sigma = 2.0, 
               theta = [1.0, 1.0],
               seed = convert(Int64, round(rand()*1e8))
               )
    # generate data
    dat = simdat_crossed(subj_n, item_n, subj_btwn = subj_btwn, both_win = both_win )

    # set contrasts
    contrasts = Dict(:age => HelmertCoding(), 
                     :condition => HelmertCoding());

    # set up model
    f = @formula dv ~ 1 + age*condition + (1|item) + (1|subj);
    m = fit(MixedModel, f, dat, contrasts=contrasts)

    # run simulation
    rng = MersenneTwister(seed);

    simulate_waldtests(
        rng, nsims, m, 
        β = beta, 
        σ = sigma, 
        θ = theta, 
        use_threads = true
    );
end
~~~~~~~~~~~~~


~~~~
mysim (generic function with 6 methods)
~~~~





Run simulations over a range of values for any parameter.

~~~~{.julia}
# varying
subj_ns = [20, 30, 40]
item_ns = [10, 20, 30]

# fixed
nsims = 1000
new_beta = [0, 0.4, 0.1, 0]
new_sigma = 2.0
new_theta = [1.0, 1.0]

d = DataFrame()

for subj_n in subj_ns
    for item_n in item_ns
        s = mysim(subj_n, item_n, nsims, new_beta, new_sigma, new_theta);
        pt = power_table(s)
        pt[!, :item_n] .= item_n
        pt[!, :sub_n] .= subj_n
        append!(d, pt)
    end
end

# save the data in long format
CSV.write("sim/power.csv", d)

# spread the table for easier viewing
unstack(d, :coefname, :power)
~~~~~~~~~~~~~


~~~~
9×6 DataFrames.DataFrame. Omitted printing of 1 columns
│ Row │ item_n │ sub_n │ (Intercept) │ age: Y   │ age: Y & condition: B │
│     │ Int64  │ Int64 │ Float64⍰    │ Float64⍰ │ Float64⍰              │
├─────┼────────┼───────┼─────────────┼──────────┼───────────────────────┤
│ 1   │ 10     │ 20    │ 0.07        │ 0.257    │ 0.059                 │
│ 2   │ 10     │ 30    │ 0.081       │ 0.345    │ 0.055                 │
│ 3   │ 10     │ 40    │ 0.088       │ 0.428    │ 0.046                 │
│ 4   │ 20     │ 20    │ 0.067       │ 0.274    │ 0.049                 │
│ 5   │ 20     │ 30    │ 0.066       │ 0.338    │ 0.064                 │
│ 6   │ 20     │ 40    │ 0.054       │ 0.413    │ 0.048                 │
│ 7   │ 30     │ 20    │ 0.045       │ 0.289    │ 0.048                 │
│ 8   │ 30     │ 30    │ 0.058       │ 0.344    │ 0.055                 │
│ 9   │ 30     │ 40    │ 0.061       │ 0.438    │ 0.058                 │
~~~~





## Convert this file 

~~~~{.julia}

# using Weave

# convert to html
# weave("simulation_tutorial.jmd")

# convert to a python notebook
# convert_doc("simulation_tutorial.jmd", "simulation_tutorial.ipynb")

# convert to md for README
# weave("simulation_tutorial.jmd", doctype="pandoc", out_path = "README.md")
~~~~~~~~~~~~~

