# MathTB2D
A Mathematica Paclet for constructing Tight-binding Hamiltonian in 2D


## How it works?

The essential part of this Paclet is constructing the Tight-binding Hamiltonian in the reciprocal space for **any** 2D lattice structure. It also allows for calculating the band spectrum. The input includes Bravis vectors for the lattice and the specific orbital informations (sites and numbers). After defining the lattice structure, the user may customize the hopping between different site and the onsite enenrgies.

## Installation

Download and put the ```MathTB2D.wl``` file in the following path
```
Mathematica\Contents\AddOns\ExtracPackages
```
or somewhere your Mathematica could find ($Path).

To load the package:
```
Needs["MathTB2D`"]
```
or one use following command, which will reload the package (overwite the variables) every time being called

```Get["MathTB2D`"] ```

The above command is equivalent with
```
<<"MathTB2D`"
```
After the package is sucessfully loaded, type ```?MathTB2D`* ``` to see all variables and available functions.

**Note:** We reserve the variables $kx$ and $ky$ for the crystal momentum. It's recommended one should not modify these two variables and avoid duplicate names for variables.

## Functions

 * ``` SetModel[lattice_, sitePositions_, orb_]```
This function sets up the model for which the tigh-binding model will be applied. The valid input ```lattice``` should contain two vectors (lists) as ```lattice={v1,v2}``` with ```v1={x1,y1}``` and ```v2={x2,y2}``` being the noncollinear 2D Bravis vectors specifying repeated unit cell of the lattice. The input ```sitePositions={s1,s2,...}``` should contain a list for the site informations for each orbital within the unit cell. The i-th orbital site is specified by relative coordinates: ```s1={s1[[1]], s1[[2]]}``` as ```s1[[1]]*v2+ s1[[2]]*v2```. Thus, s1[[1]] and s1[[2]] must satisfy ``` 0<=s1[[1]], s1[[2]] <1```. The last input ```orb={o1,o2,o3}``` contains the number of orbitals at each site with ```o1,o2,o3``` being the positive integers. Note that the length of ```orb``` must match the length of ```sitePositions```.

 **Note:** The "orbitals" does not have to represent the actual atomic orbitals. They could represent any physical degree of freedom associated with the corresponding site. For example, ```o1=o2=2``` could represent that we have two sites in one unit cell and each site has two orbitals. But they could also represent that we only have one orbital at each site but with additional spin up and spin down degree of freedom (spinful calculation). The physical meaning of "orbital" completely depends on the user's problem. Designed in this way, the spin-orbital coupled could be conveniently added by the "orbital-resolved" hoopping.
 <br/><br/>
 

 * ``` GetModel[]```

This function returns the model parameters under the current setup.
<br/><br/>


 * ``` CheckLatticeStructure[]```

This function retunrs ```True``` is the defined model has a valid lattice structure. Othewise, return ```False```.
<br/><br/>


 * ``` PlotLatticeStructure[]```

This function returns a "Manipulate plot" for visualizaing the 2D lattice structure. There are three adjustable parameters in the plot, ```n1,n2,NN``. The first determines how many repeated unit cells are plotted and the last one ```NN``` controls the demonstration of bond. For ```NN=0```, only sites are plotted (as dots in different color). For ```NN=1```, only bands that connect the nearest-neighbor sites are plotted. For ```NN=2```, bands that connect both the nearest-neighbor sites and next-nearest-neighbor sites are plotted. One can increase ```NN``` up to ```NN=5``` while the bonds are plotted in black color with decreasing weight.
<br/><br/>



 * ``` AddHoppings[t_, i_, j_, R_]```

This function adds the hopping terms for the defined lattice. The input ```t``` is the amplitude for the hoppings from ```i-th``` site of the source unit cell to the ```j-th``` site of the target unit cell.

**Note that**: 

(1) The targe unit cell has a relative shift given by ```R``` with respect to the source unit cell. The shift given by ```R={r1,r2}``` is ```r1*v1+r2*v2```. Thus, ```r1``` and ```r2``` must be integers. The shift ```R``` could be zero, which represents the hopping is within one unit cell. Howver, when ```R={0,0}```, we must enforce that $i\neq j$, otherwise this would mean the energy at the same site. For setting the onsite energy, please use ``` AddOnsite[onsite_]``` instead (see following). 


(2) Once the hopping from the source unit cell to the targer unit cell is set, there is no need to set the hopping from the target unit cell to the source unit cell. This function will automatically do that to ensure the Hermicity. 


(3) For case of multiple orbitals at each site, the hopping amplitude ```t``` is a matrix. For example, when specifying the hopping between first and second site, ```t``` is an o1-by-o2 matrix. The spin-orbti coupling can be include by adding a complex hopping ```t```. For the single-orbital case, the hopping amplitude is a scalare but must be written in the 1-by-1 matrix form as ```t={{t0}}```. 

(4) The hopping  ```t``` matrix does not need to the a Hermitian matrix. In fact, it could be even non-square for the hoppings between sites with different number of orbitals. Howver, the final Hamiltonian will always be Hermitian.

with additional 
<br/><br/>



 * ``` AddOnsite[onsite_]```

This function adds onsite energy to the model. The input ```onsite``` must be a list with real number and has the same length as the sum of all the orbital numbers at each sites within the unit cell (```Ttoal[orb]```). In another word,  The input ```onsite``` must has the same length of the Hamiltonian and it's added directly in the diagonal elements of the Hamiltonian. The first ```o1``` elements of ```onsite``` list specify the onsite energies of the corresponding ```o1``` orbitals at the first site of the unit cell. The second ```o2``` elements of ```onsite``` list specify the onsite energies of the corresponding ```o2``` orbitals at the second site of the unit cell. And so on so forth.
<br/><br/>


 * ``` GetHK[]```

This function returns the tight-binding Hamiltonian in the reciprocal space under the current setup of hoppings and onsite energies for the desgined lattice.
<br/><br/>


 * ``` ResetOnsite[]```

This function clears all the onsite energy.
<br/><br/>


 * ``` ResetHoppings[]```

This function clears all the hopping terms.
<br/><br/>


 * ``` PlotBands[xrange_, yrange_, Nx_, Ny_] ```

This function plots the band spectrum of the designed tight-binding model. The inputs ```xrange={kxMin, kxMax}``` and ```yrange={kyMin, kyMax}``` contain the limit of crystal momentum $kx$ and $ky$ in the plotting. The inputs ```Nx``` and ```Ny``` control the sampling points along the kx and ky direction, respectively. Thus, one will have a $Nx\times Ny$ square mesh grid for plotting the band spectrum. 
<br/><br/>



## Examples
In the following, I demonsate how my package can be applied to honeycomb lattice and the Kagome lattice. In the first case, we can see that in the band structure the Dirac cones around the K and K' valleys are successfully retrieved. In the second case, we see an additional flat bands. The notebook for these two examples can be found in the Demo folder. 






## Updates

* Version-1.0 2025/07/31
  
  First version of MathTB2D!
