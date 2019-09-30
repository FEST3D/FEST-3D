title: Test

Since it is not possible to perform unit tests with the current status of FEST-3D,
we have defined a few integrated test cases. Using this integrated method, 
you can test the implementation of various flux-schemes, higher-order methods, and turbulence models.
Any new scheme added to FEST-3D can also be tested using the same framework.
A python script, `Test.py`, is provided in the `test/` directory to run the integrated tests
with a particular flux scheme, higher-order method, and turbulence model.
Use the following command to run all tests in the `test/` directory.

### How to run Test.py script
```
$cd <FEST-3D/root/directory>/test/
$python Test.py <arg1> <arg2> <arg3>
```

- <arg1\> : Flux scheme
    - allowed options:
        - ausm
        - ldfss0
        - ausmUP
        - ausmP
        - slau

- <arg2\> : Higher order method
    - allowed options:
        - muscl
        - ppm
        - weno

- <arg3\> : Turbulence model
    - allowed options:
        - sst
        - sst2003
        - kkl
        - sa

Examples:

1. $python Test.py ausm   muscl sst
2. $python Test.py slau   weno  sa
3. $python Test.py ausmUP ppm   kkl


You will see following output on the screen:
```

  ----- Integrated Tests Started -----  
Total two processes will be used with MPICH library
Running Test number 1  --->  Subsonic flow over a smooth bump
Running Test number 2  --->  Laminar flow over a flat plate
Running Test number 3  --->  Turbulent flow over a flat plate
 ----- All tests completed -----
 
Tests passed:  3 out of 3
Check test summary in 'Report.txt' file.

```

The integrated tests use two processes with MPICH library. Three different test cases are defined:

1. Inviscid test case: [Subsonic flow over a 2D smooth bump](./05_tutorials/02_2dbump.html).
2. Laminar test csae: [Laminar flow over a flat plate](./05_tutorials/04_LamFp.html).
3. Turbulent test case: [Fully turbulent flow over a flat plate](./05_tutorials/05_TurbFp.html).

| S.No.&nbsp;&nbsp;&nbsp; 	| Test case &nbsp;&nbsp;&nbsp;	| Expected value &nbsp;&nbsp;&nbsp; | Tolereance % 	|
|:-------	:|:-----------	|:------------------------------	|:------------	:|
| 1     	| Inviscid  	| Change in entropy = 0.0%      	| 0.1%       	|
| 2.    	| Laminar   	| Cofficient of drag = 0.00133 	| 1%         	|
| 3.    	| Turbulent 	| Cofficient of drag = 0.00290 	| 2%         	|


The percentage change or error(Tolerance) in entropy is calculated using following expression:
<figure>
  <div style="display: flex; flex-wrap: wrap; justify-content: center; align-items:center">
    <img src="|media|/error.png" alt="Domain" style="width:250px">
  </div>
</figure>

Here, **S** is entropy, and **V** is volume. Subscript _cell_ is used for cell variable,
 _infinity_ is used for freestream quantity, and _total_ is used for to whole domain value.


For more details about domain, boundary conditions and flow conditions of these test cases, 
check the Tutorial section. The solver setup (domain, grid, flow and boundary conditions) is same as listed in the separate tutorials.

Once the tests are complete, you can check the test summary in `Report.txt` file.
``` 
Ran Test number 1  --->  Subsonic flow over a smooth bump
 __________Report__________  
 ---------- Inviscid Test case: Smooth Bump ---------- 
 Flux Scheme        : ausm
 Higher order method: muscl
 Turbulence model   : none
 Expected Change in entropy           : 0.000E+00
 Calculated relative change in entropy: 1.046E-06
 Difference                           : 1.046E-04 %
 Allowed Tolerance                    : 0.1 %
------------ >>> Test Passed  <<< --------------


Ran Test number 2  --->  Laminar flow over a flat plate
  __________Report__________   
 ------ Laminar Test case: Flat plate ------ 
 Flux Scheme        : ausm
 Higher order method: muscl
 Turbulence model   : none
 Expected drag coeffcient    : 1.330E-03
 Calculated drag coefficient : 1.329E-03
 Difference                  : 4.638E-02 %
 Allowed Tolerance           : 1 %
------------ >>> Test Passed  <<< --------------


Ran Test number 3  --->  Turbulent flow over a flat plate
  __________Report__________   
 ------ Turbulent Test case: Flat plate ------ 
 Flux Scheme        : ausm
 Higher order method: muscl
 Turbulence model   : sst
 Expected drag coeffcient    : 2.900E-03
 Calculated drag coefficient : 2.873E-03
 Difference                  : 9.312E-01 %
 Allowed Tolerance           : 2 %
------------ >>> Test Passed  <<< --------------
```

