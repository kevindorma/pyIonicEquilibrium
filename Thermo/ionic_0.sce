// Kev's ionic solver
// October 2017



// Algorithm

//obtain constraint info
//constraint(1:constraint) = [:,3]
//   :,1 is i, :,2 is j, :,3 is value

//define column position for each of the i and j positions
//k(i,j) = (i-1) + (j-1)*nj + 1

//initialize guesses and matrix A x = b
//initialize x_prev and x
//Get the constrains (known values).  A_known x = b

//Construct block matrix for linear equations.
//Construct block matrix for ln equations
//Construct block matrix for x -> ln x definitions  [careful, this will fail for initial conditions of zero]
//Construct block matrix for constraints

//while not converged

//Update the block matrix for x -> ln x definitions
//Assemble matrix equation A x = b

//Solve A x = b
//if x is not close to x_prev then set x_prev = x, set converged = FALSE
//otherwise set converged = TRUE
//end while

// set reaction temperature
rxnT0 = 298.15
rxnT1_C = 120
rxnT1 = rxnT1_C + 273.15


// create list of species
species_w_water = ["OH-","H2S","HS-","CO2","HCO3-","CO32-","H+","S2-","H2O"];
nspecies_w_water = max(size(species_w_water))
species = species_w_water(1:nspecies_w_water-1)
nspecies = max(size(species))

// now fill in the properties of each species
Hform = zeros(nspecies_w_water,1)    // std state 298 K
Cpmolar = zeros(nspecies_w_water,1) // this is at std state of 298 K
// NIST webprop is a good place to find this information
Hform0 = [-230;	-20.6;	-18;	-393;	-692;	-677;	0;	33;	-286]


// I might not like the convention used in HnH.  Is 
// the basis is Cp for H+ is zero.
// the absolute value for H+ is near -71 J/mol.K (Ion Properties, page 119)
//  This is a great resource for ion heat capacities
// Table 9 has OH- as -69
// S2- = -184
// CO32- = -159
// HS-  = -23
// HCO3-  = 18


for i = 1:nspecies_w_water
    select species_w_water(i)
    case "OH-"
        Cpmolar(i) = -69 
    case "H2S"
        Cpmolar(i) = 178 
    case "HS-"
        Cpmolar(i) = -23
    case "CO2"
        Cpmolar(i) = 243
    case "HCO3-"
        Cpmolar(i) = 18
    case "CO32-"
        Cpmolar(i) = -159
    case "H+"
        Cpmolar(i) = -71
    case "S2-"
        Cpmolar(i) = -184 
    case "H2O"
        Cpmolar(i) = 4.2 * 18
    end;
    
end





// create list of reaction names
rxname = ["bisulfide";"bicarbonate";"carbonate";"water";"disulfide"]
nrxn = max(size(rxname))

//H2S + NaOH = NaHS + H2O
//CO2 + NaOH = NaHCO3
//NaHCO3 + NaOH = Na2CO3 + H20
//H+ + OH- = H2O
//NaHS + NaOH = Na2S + H2O

//now create a stoichiometry matrix for the reactions
// I can make this more elegant, but that is for later.
// define Hformation(nspecies,1), DHrxn_std(nrxn,1), molarCp(nspecies,1), DHrxn_t(nrxn,1), and revise the eqm constants.
// just hack this, no fancy stuff
rxn_w_water = [1, 1, -1, 0, 0, 0, 0, 0, -1;
        1, 0, 0, 1, -1, 0, 0, 0, 0;
        1, 0, 0, 0, 1, -1, 0, 0, -1;
        1, 0, 0, 0, 0, 0, 1, 0, -1.
        1, 0, 1, 0, 0, 0, 0, -1, -1]
rxn = rxn_w_water(:,1:nspecies)
// Keqm values, vector
Keqm = [1.85E+08;
    1.62E+09;
    4.20E+04;
    1.93E+15;
    1.30E+00]
lnKeqm = log(Keqm)

// and do this at std state
Keqm0 = [1.10E+07;
    4.20E+07;
    4.80E+03;
    1.00E+14;
    1.00E+00];

lnKeqm0 = log(Keqm0)


// and now create the standard heat of reaction at 298 K
Hrxn0 = rxn_w_water * Hform0 * 1000
// and now extrapolate to a higher temperature
// this requires estimating the heat of reaction at the midpoint from 298 to desired temperature
// couple of ways to do this.
// use the std heat of reaction, and add Delta Cp * 0.5 * (T1 - T0)
Hrxn1 = Hrxn0 + rxn_w_water * Cpmolar * 0.5 * (rxnT1 - rxnT0)

// and now we need to extrapolate the K to the new temperature
// do this later.


// initial concentration, exclude water
init_conc = [0.00E+00, 3.04E-03, 0.00E+00, 1.02E-04, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00]'

soln = ones(2*nspecies,1)
conc = ones(nspecies,1)
lnconc = ones(nspecies,1)
prev_conc = ones(nspecies,1)
prev_lnconc = ones(nspecies,1)


for counter = 1:30
    
// start assembling tghe matrices
// balance equations
leftbal = eye(nspecies,nspecies)
midbal = rxn'
rightbal = zeros(leftbal)
b_bal = init_conc

//linearized form of the eqm relations for each reaction

leftlineqm = zeros(rxn)
midlineqm = zeros(nrxn,nrxn)
rightlineqm = -rxn;
b_lineqm = lnKeqm;

// linearize the log species
leftlinlog = eye(nspecies,nspecies)
midlinlog = zeros(rxn')
rightlinlog = diag(-exp(prev_lnconc))
b_linlog = exp(prev_lnconc) .* (ones(prev_lnconc)-prev_lnconc)

// and assemble
Awhole =  [ leftbal, midbal, rightbal;
            leftlineqm, midlineqm, rightlineqm;
            leftlinlog, midlinlog, rightlinlog]
bwhole = [b_bal ; b_lineqm ; b_linlog];

soln = Awhole\bwhole
conc = soln(1:nspecies)
extent = soln(nspecies+1:nspecies+nrxn)
lnconc = soln(nspecies+nrxn+1:2*nspecies+nrxn)

mysoln(counter,:)=soln'

residual(counter) = norm(lnconc-prev_lnconc)
old_soln(counter,1:nspecies)=lnconc'

prev_soln = soln
prev_conc = conc
prev_lnconc = lnconc

// and it works
// yay

end;
