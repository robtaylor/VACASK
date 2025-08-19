![Xschem showing the rfosc.sch example.](rfosc.png)

# Using Xschem with VACASK

Download the latest Xschem sources from [https://github.com/StefanSchippers/xschem](https://github.com/StefanSchippers/xschem). Most packages available in various Linux distributions lack VACASK support because they are too old. Build and install Xschem (by default it will be installed in /usr/local). 
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
|cccs, ccvs                  |quote `vnam` value, escape quotes (e.g. `vnam=\"vctl\"`) |
|res, capa, capa-2, ind      | |
|k                           |quote `L1` and `L2` value, escape quotes |
|npn, pnp, njfet, pjfet      | |
|nmos, pmos                  | |
|nmos3, pmos3                |uses `m` instead of `$mfactor` (assumes device is a subcircuit)|
|nmos4, pmos4                | |
|nmos-sub, pmos-sub          | |
|nmos4_depl                  | |
|pmoshv4, pmosnat            |uses `number` instead of `$mfactor` (assumes device is a subcircuit)|
|code, code_shown            | |
|netlist, netlist_not_shown  | |
|netlist_at_end              | |
|netlist_not_shown_at_end    | |
|simulator_commands          |set `simulator` to `vacask` |
|simulator_commands_shown    |set `simulator` to `vacask` |
|param                       | |

To make a symbol VACASK compatible set the `spectre_format` attribute in the .sym file. The `spectre_device_model` attribute specifies the VACASK netlist commands that will be included in the netlist once for that particular type of element (i.e. `load` and `model` commands for fundamental elements, like independent sources, controlled sources, and passives). Before creating a netlist, select the `spectre` netlist mode in the `Options/Netlist format` menu. 

An Xschem configuration file ([xschemrc](xschemrc)) is provided that sets the correct netlist type, configures the directory for storing the netlists and running simulations, and sets up the external editor. VACASK must be in the system path if you want to run simulations from Xschem. 

VACASK rawfiles always contain a single plot. When loading raw files with the `xschem raw_read <rawfile.raw> <analysis type>` command do not specify the analysis type. This way Xschem will load the first plot from the rawfile and will not try to match the `Plotname:` field in the file. The same can be achieved by selecting the `Waves/Load first analysis found` option from the main menu, as well as, by clicking the `View results` icon. 

A simple example is provided ([rfosc.sch](rfosc.sch)). After you open it, select `Netlist` from the main menu, followed by `Simulate`. When the simulation is finished select `Waves/Load first analysis found` in the main menu and choose file `tran1.raw`. The circuit's response will be plotted in the embedded graph. 

