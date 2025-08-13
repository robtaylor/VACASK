# Converting IHP Open PDK for use with VACASK

Clone IHP SG13G2 PDK
```
git clone https://github.com/IHP-GmbH/IHP-Open-PDK
```
Now you have a directory named `IHP-Open-PDK`. 

Set environmental variables
* `PDK_ROOT` to the directory where you downloaded the PDK (e.g. `/home/myname/IHP-Open-PDK`) and
* `PDK` to `ihp-sg13g2`. 

You will need the path to VACASK's Python scripts. If you don't know where these scripts are, type
```
vacask -dp
```
and look for "Python path addition". Suppose the python path addition is `/usr/local/lib/vacask/python`. Type
```
PYTHONPATH=/usr/local/lib/vacask/python python3 -m sg13g2tovc
```

The converter will process the Ngspice models and 
* create directory `ihp-sg13g2/libs.tech/vacask/models` with the converted models
* create directory `ihp-sg13g2/libs.ref/sg13g2_stdcell/vacask` with the converted standard cells
* create directory `ihp-sg13g2/libs.ref/sg13g2_io/vacask` with the converted I/O cells
* create a VACASK config file `ihp-sg13g2/libs.tech/vacask/.vacaskrc.toml`
  (copy this file to the directory where your top level netlist is located)
* patch xschem symbols with VACASK netlisting pattern (currently only the sg13g2_pr directory is processed)
* compile the Verilog-A models that are provided with the PDK and place the 
  resulting .osdi files in `ihp-sg13g2/libs.tech/vacask/osdi`
* create an include file `ihp-sg13g2/libs.tech/vacask/models/sg13g2_vacask_common.lib`
  (always include this file along with all other PDK files)

Some examples and a [.vacaskrc.toml](.vacaskrc.toml) configuration file are available in this directory. Just set the `PDK_ROOT` and the `PDK` environmental variables and run examples with VACASK. 

Xschem support is under development. Use the latest development version available on [Github](https://github.com/StefanSchippers/xschem) (you will have to compile it yourrself). Copy file `ihp-sg13g2/libs.tech/xschem/xschemrc` to the directory where you want to start xschem. It will load the IHP PDK's extensions and symbols at startup. Select `Options/Netlist format/Spectre` before you dump the netlist or run simulations. 