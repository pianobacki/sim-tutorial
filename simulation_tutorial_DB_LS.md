Simulation Tutorial
===================

Based on Lisa DeBruine's work.

Updated by Lisa Schwetlick & Daniel Backhaus


### Load the packages we'll be using in Julia
```@example Main
using MixedModels        # run mixed models
using MixedModelsSim     # simulation functions for mixed models
using RCall              # call R functions from inside Julia
using DataFrames, Tables # work with data tables
using StableRNGs         # random number generator
using CSV                # write CSV files
using Markdown
using Statistics         # basic math funcions
using LinearAlgebra      # not used yet, for specifing θ
```

### Define number of iterations

Set to a low number for test, high for production
```@example Main
nsims = 1000 
```
<br/><br/>
# Simulate from existing data

Load existing data:
```@example Main
kb07 = MixedModels.dataset("kb07");
```

Set contrasts:
```@example Main
contrasts = Dict(:spkr => HelmertCoding(), 
                 :prec => HelmertCoding(), 
                 :load => HelmertCoding());
```

Define formula:
```@example Main
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );
```

Fit the model
```@example Main
kb07_m = fit(MixedModel, kb07_f, kb07, contrasts=contrasts);
print(kb07_m)
```

## **Simulate from existing data with same parameters**
 
Use the `parameparametricbootstrap()` function to run `nsims` iterations of data sampled using the parameters from `kb07_m`. 
Set up a random seed to make the simulation reproducible. You can use your favourite number.

To use multithreading, you need to set the number of cores you want to use. 
E.g. in Visual Studio Code, open the settings (gear icon in the lower left corner or cmd-,) and search for "thread". 
Set `julia.NumThreads` to the number of cores you want to use (at least 1 less than your total number).


Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Run nsims iterations:
```@example Main
kb07_sim = parametricbootstrap(rng, nsims, kb07_m, use_threads = false);
```

**Try**: Run the code above with or without `use_threads = true`.

Convert p-values to dataframe and save it as CSV
```@example Main
kb07_sim_df = DataFrame(kb07_sim.coefpvalues);
CSV.write("kb07_sim.csv", kb07_sim_df);
```

Have a look:
```@example Main
print(first(kb07_sim_df, 8))
```

### Power calculation

The function `power_table()` from `MixedModelsSim` takes the output of `parametricbootstrap()` and calculates the proportion of simulations where the p-value is less than alpha for each coefficient. 
You can set the `alpha` argument to change the default value of 0.05 (justify your alpha ;).

```@example Main
ptbl = power_table(kb07_sim)
```

You can also do it manually:
```@example Main
kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:]

mean(kb07_sim_df[kb07_sim_df.coefname .== Symbol("prec: maintain"),:p] .< 0.05)
````

For nicely displaying, you can use pretty_table:
```@example Main
pretty_table(ptbl)
```


## **Simulate data with changed parameters without touching the existing data**


Let's say we want to check our power to detect effects of spkr, prec, and load 
that are half the size of our pilot data. We can set a new vector of beta values 
with the `β` argument to `parametricbootstrap()`.


Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Specify β:
```@example Main
new_beta = kb07_m.β
new_beta[2:4] = kb07_m.β[2:4]/2
```

Run nsims iterations:
```@example Main
kb07_sim_half = parametricbootstrap(rng, nsims, kb07_m, β = new_beta, use_threads = false);
```

### Power calculation

```@example Main
power_table(kb07_sim_half)
```

<br/><br/>


# Simulate data from scratch

If we simulate data from scratch, next to subject and item number, we can manipulate the arguments `β`, `σ` and `θ`.
Lets have a closer look at them, define their meaning and we will see where the corresponding values in the model output are.

### **Beta**
`β` are our effect sizes. If we look again on our LMM summary from the `kb07`-dataset `kb07_m`
we see our four `β` under fixed-effects parameters in the `Coef.`-column. 

```@example Main
kb07_m 
kb07_m.β
```

### **Sigma**
`σ` is the residual-standard deviation listed under the variance components. 
```@example Main
kb07_m 
kb07_m.σ
```

### **Theta**
`θ` is a more complex parameter. In a less complex model, with only intercepts for the random effects, 
or if we supress the correlations in the formula with `zerocorr()` then `θ` describes the relationship between 
the random effects standard deviation and the standard deviation of the residual term.
In our `kb07_m` example:
The residual standard deviation is `680.032`.
The standard deviation of our first variance component *`item - (Intercept)`* is `364.713`.
Thus our first `θ` is the relationship: variance component devided by residual standard deviation
364.713 /  680.032 =  `0.53631`

```@example Main
kb07_m.θ
```

We also can calculate the `θ` for variance component *`subj - (Intercept)`*. 
The residual standard deviation is `680.032`.
The standard deviation of our variance component *`subj - (Intercept)`* is `298.026`.
Thus the retated θ is the relationship: variance component devided by residual standard deviation
298.026 /  680.032 =  `0.438252`

```@example Main
kb07_m.θ
```

We can not calculate the `θ` for variance component *`item - prec: maintain`* yet, because it includes the correlation of 
*`item - prec: maintain`* and *`item - (Intercept)`*. 
The `θ` vector is the flattened version of the variance-covariance matrix - a lowertrinangular matrix.
The on-diagonal elements are just the standard deviations (the `σ`'s), If all off-diagonal elements are zero, we can use our
calculation above. The off-diagonal elements are covariances and correspond to the correlations (the `ρ`'s). 
If they are unequal to zero, as it is in our `kb07`-dataset, we cannot recreate the variance-covariance matrix having the model output.
We just take it from the model we already fitted

See the two inner values:
```@example Main
kb07_m.θ
```

# NEED HELP HERE, IS THEIR ANY WAY TO DEFINE VARIANCE-COVARIANCE MATRIX YET ?
<!---
vc = VarCorr(kb07_m)
vc.σρ
kb07_m.λ
kb07_m.λ[1]
kb07_m.λ[2]


λitem = LowerTriangular(diagm([1.3, 0.35, 0.75]))
λsubj = LowerTriangular(diagm([1.5, 0.5, 0.75]))

isapprox(kb07_m.θ,  [flatlowertri(kb07_m.λ[1]); flatlowertri(kb07_m.λ[2])])

[flatlowertri(λitem); flatlowertri(λsubj)]


# make a lower triangular matricis
re_item = create_re(0.5363168233715857,0)
re_item[4]=0.37133693708531357

re_subj = create_re(0.4382528181348316)
# make the compact form out of it = is equal to θ 
vcat( flatlowertri(re_item), flatlowertri(re_subj) )
--->

Having this knowledge about the parameter we can now **simulate data from scratch**

The `simdat_crossed()` function from `MixedModelsSim` lets you set up a data frame with a specified experimental design. 
For now, it only makes fully balanced crossed designs!, but you can generate an unbalanced design by simulating data for the largest cell and deleting extra rows. 

Firstly we will set an easy design where `subj_n` subjects per `age` group (O or Y) respond to `item_n` items in each of two `condition`s (A or B).

Your factors need to be specified separately for between-subject, between-item, and within-subject/item factors using `Dict` with the name of each factor as the keys and vectors with the names of the levels as values.
<br/><br/>
Define factors in a dict.
Put between subject factors in a dict:
```@example Main
subj_btwn = Dict("age" => ["O", "Y"])
```

There are no between-item factors in this design so you can omit it or set it to nothing. 
Note if you have factors which are between subject and between item, put them in both dicts.
```@example Main
item_btwn = nothing
```

Put within-subject/item factors in a dict:
```@example Main
both_win = Dict("condition" => ["A", "B"])
```

Define subject and item number:
 ```@example Main
subj_n = 10
item_n = 30
```

Simulate data:
```@example Main
dat = simdat_crossed(subj_n, 
                     item_n, 
                     subj_btwn = subj_btwn, 
                     item_btwn = item_btwn, 
                     both_win = both_win);
```

Have a look:
```@example Main
first(DataFrame(dat),8)
```
The values we see in the column `dv` is just random noise.

Set contrasts:
```@example Main
contrasts = Dict(:age => HelmertCoding(), 
                 :condition => HelmertCoding());
```

Define formula:
```@example Main
f1 = @formula dv ~ 1 + age * condition + (1|item) + (1|subj);
```

Fit the model:
```@example Main
m1 = fit(MixedModel, f1, dat, contrasts=contrasts)
print(m1)
```

Because the `dv` is just random noise from N(0,1), there will be basically no subject or item random variance, 
residual variance will be near 1.0, and the estimates for all effects should be small. 
Don't worry, we'll specify fixed and random effects directly in `parametricbootstrap()`. 


Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Specify `β`, `σ`, and `θ`, we just made up this parameter:
```@example Main
new_beta = [0., 0.25, 0.25, 0.]
new_sigma = 2.0
new_theta = [1.0, 1.0]
```

Run nsims iterations:
```@example Main
sim1 = parametricbootstrap(rng, nsims, m1, 
                        β = new_beta, 
                        σ = new_sigma, 
                        θ = new_theta,
                        use_threads = false);
```

### Power calculation

```@example Main
ptbl= power_table(sim1)
```

For nicely displaying it, you can use pretty_table:
```@example Main
pretty_table(ptbl)
```

<br/><br/>
# Recreate the `kb07`-dataset from scratch

To also play with the number of subjects and items in our power calculation, we need to recreate the 
design of an given experiment. Lets start with a simulated replication of the `kb07`-dataset.


Define subject and item number:
```@example Main
subj_n = 56
item_n = 32
``` 

Define factors in a dict:
```@example Main
subj_btwn = nothing
item_btwn = nothing
both_win = Dict("spkr" => ["old", "new"],
                "prec" => ["maintain", "break"],
                "load" => ["yes", "no"])
```

Simulate data:
```@example Main
fake_kb07 = simdat_crossed(subj_n, item_n, 
                     subj_btwn = subj_btwn, 
                     item_btwn = item_btwn, 
                     both_win = both_win);
```

Make a dataframe:
```@example Main
fake_kb07_df = DataFrame(fake_kb07)
```

Have a look:
```@example Main
first(fake_kb07_df,8)
```

# NEED HELP - IS THERE AN EASIER WAY ?

Change sorting for later selection:
```@example Main
fake_kb07_df = unstack(fake_kb07_df, :item, :dv);   # makes wide format  
fake_kb07_df = stack(fake_kb07_df, 5:ncol(fake_kb07_df), variable_name = :item);  # makes long format
rename!(fake_kb07_df, :value => :dv)

fake_kb07_df = unstack(fake_kb07_df, :subj, :dv)
fake_kb07_df = stack(fake_kb07_df, 5:ncol(fake_kb07_df), variable_name = :subj);
```

Rename the `dv` to its original name in the `kb07`-dataset:
```@example Main
rename!(fake_kb07_df, :value => :rt_trunc)
``` 

Write a CSV:
```@example Main
CSV.write("fake_kb07_df.csv", fake_kb07_df);
```

# NEED HELP, is that correct? or is it possible to do it in simdat_crossed?

Our original design is not fully crossed. Every subject saw an image only once, thus in one of eight possible conditions. To simulate that we only keep one of every eight lines.

Define an index which represets a random choise of one of every eight rows:
```@example Main
Z= convert(Int64,(length(fake_kb07)/8))
idx = rand(rng, 1:8 , Z)
A = repeat([8], inner=Z-1)
A = append!( [0], A )
A = cumsum(A)
idx = idx+A
```

Reduce the fully crossed design to the original experimental design:
```@example Main
fake_kb07_df= fake_kb07_df[idx, :]
```

Set contrasts:
```@example Main
contrasts = Dict(:spkr => HelmertCoding(), 
                 :prec => HelmertCoding(), 
                 :load => HelmertCoding());
``` 

Define formula, same as above:
```@example Main
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );
```

Fit the model:
```@example Main
fake_kb07_m = fit(MixedModel, kb07_f, fake_kb07_df, contrasts=contrasts);
print(fake_kb07_m)
```

Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Specify `β`, `σ`, and `θ`, we changed the parameter to the values we already know from the model of the existing dataset:

```@example Main
new_beta = [2181.85, 67.879, -333.791, 78.5904]
new_beta = kb07_m.β

new_sigma = 680.032
new_sigma = kb07_m.σ

new_theta = [0.5363168233715857,
           -0.25981993107379,
            0.2653016002105174,
            0.4382528181348316]
new_theta = kb07_m.θ
```


Run nsims iterations:
```@example Main
fake_kb07_sim = parametricbootstrap(rng, nsims, fake_kb07_m, 
                        β = new_beta, 
                        σ = new_sigma, 
                        θ = new_theta,
                        use_threads = false);
```

### Power calculation

```@example Main
power_table(fake_kb07_sim)
```

Compare to the powertable from the existing data:
```@example Main
power_table(kb07_sim)
```

We have successfully recreate the powersimulation of an existing dataset from scratch. This has the 
advantage, that we now can iterate over different number of subjects and items.

<br/> <br/>

# Loop over sample sizes

Run simulations over a range of values for any parameter.
Therefore we first define every fixed things outside the loop.

Define factors in a dict:
```@example Main
subj_btwn = nothing
item_btwn = nothing
both_win = Dict("spkr" => ["old", "new"],
                "prec" => ["maintain", "break"],
                "load" => ["yes", "no"]);
```

Set contrasts:
```@example Main
contrasts = Dict(:spkr => HelmertCoding(), 
                 :prec => HelmertCoding(), 
                 :load => HelmertCoding());
``` 

Set random seed for reproducibility:
```@example Main
rng = StableRNG(42);
```

Define formula, same as above:
```@example Main
kb07_f = @formula( rt_trunc ~ 1 + spkr+prec+load + (1|subj) + (1+prec|item) );
```

Specify `β`, `σ`, and `θ`, we changed the parameter to the values we already know from the model of the existing dataset:

```@example Main
new_beta = [2181.85, 67.879, -333.791, 78.5904]
new_beta = kb07_m.β

new_sigma = 680.032
new_sigma = kb07_m.σ

new_theta = [0.5363168233715857,
           -0.25981993107379,
            0.2653016002105174,
            0.4382528181348316]
new_theta = kb07_m.θ
```



Define subject and item numbers as arrays:
```@example Main
sub_ns = [20, 30, 40];
item_ns = [10, 20, 30];
``` 

Make an emty dataframe:
```@example Main
d = DataFrame();
```

Run the loop:
```@example Main
for sub_n in sub_ns
    for item_n in item_ns

    # Make fully crossed data:
    fake_kb07 = simdat_crossed(subj_n, item_n, 
                     subj_btwn = subj_btwn, 
                     item_btwn = item_btwn, 
                     both_win = both_win);


    # Reduce the fully crossed design to the original experimental design:
    fake_kb07_df = DataFrame(fake_kb07);
    fake_kb07_df = unstack(fake_kb07_df, :item, :dv);
    fake_kb07_df = stack(fake_kb07_df, 5:ncol(fake_kb07_df), variable_name = :item);

    rename!(fake_kb07_df, :value => :dv);
    fake_kb07_df = unstack(fake_kb07_df, :subj, :dv)
    fake_kb07_df = stack(fake_kb07_df, 5:ncol(fake_kb07_df), variable_name = :subj)
    rename!(fake_kb07_df, :value => :rt_trunc)

    Z= convert(Int64,(length(fake_kb07)/8))
    idx = rand(rng, 1:8 , Z)
    A = repeat([8], inner=Z-1)
    A = append!( [0], A )
    A = cumsum(A)
    idx = idx+A

    fake_kb07_df= fake_kb07_df[idx, :]
    

    # Fit the model:
    fake_kb07_m = fit(MixedModel, kb07_f, fake_kb07_df, contrasts=contrasts);

    # Run nsims iterations:
    fake_kb07_sim = parametricbootstrap(rng, nsims, fake_kb07_m, 
                        β = new_beta, 
                        σ = new_sigma, 
                        θ = new_theta,
                        use_threads = false);
    

    # Power calculation
    ptbl = power_table(fake_kb07_sim)
    ptdf = DataFrame(ptbl)
    ptdf[!, :item_n] .= item_n
    ptdf[!, :sub_n] .= sub_n
    append!(d, ptdf)

    end
end
```

Save the powertable as CSV
```@example Main
CSV.write("power.csv", d)
```








not used: 
<!---
```@example Main
R"""
require(tidyverse, quietly = TRUE)   # for data wrangling and visualisation
""";
```

# 
# It's useful to be able to weave your file quickly while you're debugging, 
# so set the number of simulations to a relatively low number while you're 
# setting up your script and change it to a larger number when everything
# is debugged.

# toDo color and facet didnt work
```{julia;label=pleasecorrect}
R"""
d <- $kb07_sim_df

ggplot(d, aes(x=!!as.name("\u03b2"),  color= coefname))+
  geom_density(show.legend = FALSE)+
  facet_wrap(~coefname, scales = "free")

"""

# or as function 
# #### Define: ggplot_betas
# 
# This function plots the beta values returned from `parametricbootstrap` using ggplot in R.
# If you set a figname, it will save the plot to the specified file.


function ggplot_betas(sim, figname = 0, width = 7, height = 5) 

    df = DataFrame(sim.coefpvalues)

    R"""
    
    p <- $df %>%
     rename(beta = !!as.name("\u03b2"))
     ggplot(aes(x=!!as.name("\u03b2"),  color=coefname ))+
      geom_density(show.legend = FALSE)+
      facet_wrap(~coefname, scales = "free")

        if (is.character($figname)) {
            ggsave($figname, p, width = $width, height = $height)
        }

        p
    """
end



# Plot betas in ggplot. In the code editor or Jupyter notebooks, you can omit the file name to just display the figure in an external window.



# just display the image
# ggplot_betas(kb07_sim) 
# save the image to a file and display (display doesn't work in weave)
ggplot_betas(kb07_sim, "kb07_betas.png")



# In documents you want to weave, save the image to a file and use markdown to display the file. Add a semicolon to the end of the function to suppress creating the images in new windows during weaving.

# ![](fig/kb07_betas.png)
```

--->