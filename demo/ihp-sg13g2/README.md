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
* patch the Xschem configuration file `ihp-sg13g2/libs.tech/xschem/xschemrc`
* add a VACASK customization file for Xschem (`ihp-sg13g2/libs.tech/xschem/xschem-vacask`)
* compile the Verilog-A models that are provided with the PDK and place the 
  resulting .osdi files in `ihp-sg13g2/libs.tech/vacask/osdi`
* create an include file `ihp-sg13g2/libs.tech/vacask/models/sg13g2_vacask_common.lib`
  (always include this file along with all other PDK files)

Some examples and a [.vacaskrc.toml](.vacaskrc.toml) configuration file are available in this directory. Just set the `PDK_ROOT` and the `PDK` environmental variables and run examples with VACASK. 

# Xschem support

Download the latest Xschem sources from [https://github.com/StefanSchippers/xschem](https://github.com/StefanSchippers/xschem). Build and install Xschem (by default it will be installed in /usr/local). 
```
git clone https://github.com/StefanSchippers/xschem
cd xschem
./configure
make
sudo make install
```

The following symbols from the devices library support VACASK: 
|Symbols    |Comment         |
|-----------|----------------|
|vsource, isource, sqwsource | |
|diode, zener                | |
|vcvs, vccs                  | |
|cccs, ccvs                  |quote vnam value, escape quotes (e.g. `vnam=\"vctl\"`) |
|res, capa, capa-2, ind      | |
|k                           |quote `L1` and vL2` value, escape quotes |
|npn, pnp, njfet, pjfet      | |
|nmos, pmos                  | |
|nmos3, pmos3                |uses `m` instead of `$mfactor` (assumed device is a subcircuit)|
|nmos4, pmos4                | |
|nmos-sub, pmos-sub          | |
|nmos4_depl                  | |
|pmoshv4, pmosnat            |uses `number` instead of `$mfactor` (assumed device is a subcircuit)|
|code, code_shown            | |
|netlist, netlist_not_shown  | |
|netlist_at_end              | |
|netlist_not_shown_at_end    | |
|simulator_commands          |set `simulator` to `vacask` |
|simulator_commands_shown    |set `simulator` to `vacask` |
|param                       | |

All symbols provided by the PDK in the `ihp-sg13g2/libs.tech/xschem/sg13g2_pr` directory are converted for use with VACASK. 

Before you start Xschem make sure the `PDK_ROOT` and the `PDK` environmental variables are set. Set the `XSCHEM_NETLIST_TYPE` environmental variable to `spectre`. and copy the patched Xschem configuration file `ihp-sg13g2/libs.tech/xschem/xschemrc` to the directory where you intend to start Xschem. This file will load the IHP PDK's extensions and symbols at startup. Now you can start Xschem. 

A simple example is provided in file [dc_lv_nmos.sch](dc_lv_nmos.sch). 

Currently HBT operating point backannotation does not work because VACASK uses the Verilog-A version of the VBIC model which does not expose operating point data. Once the VBIC model is converted from Ngspice to Verilog-A this feature will be added. 
