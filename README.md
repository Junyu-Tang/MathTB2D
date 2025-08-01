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



 * ``` AddOnsite[onsite_]```

This function adds onsite energy to the model. The input ```onsite``` must be a list with real number and has the same length as the sum of all the orbital numbers at each sites within the unit cell (```Ttoal[orb]```). In another word,  The input ```onsite``` must has the same length of the Hamiltonian and it's added directly in the diagonal elements of the Hamiltonian.
<br/><br/>


 * ``` ResetOnsite[]```

This function clears all the onsite energy.
<br/><br/>


 * ``` ResetHoppings[]```

This function clears all the hopping terms.
<br/><br/>


 * ``` PlotBands[xrange_,yrange_,Nx_,Ny_] ```

This function clears all the hopping terms.
<br/><br/>






