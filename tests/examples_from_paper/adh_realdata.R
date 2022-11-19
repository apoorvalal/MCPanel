## Loading Source files
libreq(MCPanel)
df = fread("tests/examples_from_paper/california/california_prop99.csv")
setup = panelMatrices(df, "State", "Year", "treated", "PacksPerCapita")
attach(setup)
# mask is for observed Y0 : 1-treatment matrix
mask = 1 - W

# lambda to compute treatment effect
ATT = \(Ytil) mean(Y[-(1:N0),-(1:T0)] - Ytil[-(1:N0),-(1:T0)])

# %% Diff in Diff
DID(Y, mask) |> ATT()
# MC-NNM
mcnnm_cv(Y, mask)  |> predict_mc()  |> ATT()
# horizontal elastic net
en_mp_rows(Y, mask, num_alpha = 40) |> ATT()
# vertical elastic net
t(en_mp_rows(t(Y), t(mask), num_alpha = 40))  |> ATT()
# ADH
adh_mp_rows(Y, mask) |> ATT()
