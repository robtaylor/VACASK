# Tolerance handling in VACASK

Tolerances on unknowns are always enforced. Tolerances on residuals are enforced only for Kirchoff current law equations of non-internal device nodes. This is due to the fact the internal nodes have only a single contribution and therefore establishing a reference value for relative tolerance (specified by the `reltol` simulator option) is not possible because the contribution approaches zero as the nonlinear solver converges. Nevertheless the solver can still report failed convergence due to a small reference value as the it approaches the solution, particularly for nodes to which only a single element is connected. This can be alleviated by setting the `relrefres` parameter to `"pointglobal"` (the maximum contribution across all residuals with the same nature at the current timepoint is used as the reference value) or `"global"` (the maximum contribution across all residuals with the same nature over the whole past is used as the reference value). Another possibility is to set the `relref` option to `"allglobal"` which is equivalent to setting the unknown (`relrefsol`), residual (`relrefres`), and LTE (`relreflte`) reference value to `"pointglobal"` (if we assume that these there optiomns are set to `relref` which means that the `relref` option decides on their value). If you still experience convergence problems try turning off residual tolerance checks by setting `nr_residualcheck` to 0. 

All absolute tolerances can be scaled simultaneously with the `tolscale` simulator option. Its default value is 1. 

# Nature and discipline import from .osdi files
Each .osdi file has its own set of natures and disciplines. They apply only to the models defined in that particular file. When computing global maxima (across the same nature) in the nonlinear solver the name of the nature is the one that is used to identify the same natures across .osdi files. You can list the natures and disciplined along with all their attributes by issuing the following command in the control block.  
```
print device_file("<file name substring>")
```
This will print the details of all files whose canonical name contains the specified substring. 

```
print device_files
```
lists the canonical names of all loaded .osdi files. 

# Default tolerance handling (SPICE mode)
This is the default. You can explicitly enable it by adding the following command to the control block. 
```
options tolmode="spice"
```
The absolute tolerances applied in this mode are 

|Node type  |Quantity     |Nature name |absolute tolerance option |
|-----------|------------ |------------|--------|
|Potential  |unknown      |.voltage    |vntol   |
|           |unknown idt  |.flux       |fluxtol |
|           |residual     |.current    |abstol  |
|           |residual idt |.charge     |chgtol  |
|Flow       |unknown      |.current    |abstol  |
|           |unknown idt  |.charge     |chgtol  |
|           |residual     |.voltage    |vntol   |
|           |residual idt |.flux       |fluxtol |

Natures whose name starts with a dot are builtin SPICE natures. Implicit equations (the corresponding unknowns and residuals) are treated as potential nodes. Idt natures are applied to reactive residual contributions. Currently they are used only in the element evaluation bypass algorithm. 

# Verilog-A tolerance handling
To enable it add 
```
options tolmode="va"
```
or
```
options tolmode="mixed"
```
to your control block. The `va` mode considers only absolute tolerances defined in .osdi files. Each node of every device applies an absolute tolerance to the circuit node to which it is connected. If two devices want to apply different absolute tolerances to the same circuit node the lower value is used along with the nature from where this tolerance was taken. VACASK performs no nature compatibility checks. 

Tolerances that are applied to circuit nodes via the nodes of a device can be listed with 
```
print device("<device name>")
```
In the output there is a section named "Absolute tolerances in Verilog-A mode". In this section the absolute tolerances of the device's nodes are listed (the unknown, its integral, the residual, and the residual's integral). The names of the natures where the unknown and the residual tolerances are taken from are listed in the "Nodes" section. The tolerances for the integral of the unknown and the residual are taken from the idt attribute of the unknown/residual nature. If the idt attribute is missing the same tolerance is used as for the unknown/residual. 

Some nodes do not apply tolerances because OpenVAF currently does not expose their natures. These nodes are
- implicit equations created from ddt() and idt() operators
- implicit equations created manually in the module
- implicit equations of switch branches (apply only a tolerance to the corresponding unknown)

Due to this some circuit unknowns/residuals may not have any tolerance assigned. For these unknowns/residuals tolerance checks are performed. 

Builtin devices (independent and controlled linear voltage/current sources) apply the standard SPICE tolerances (see .voltage/.current/.flux/.charge above). 

When `tolmode="mixed"` is set a mixture of Verilog-A and SPICE tolerances is used. Unknowns/residuals that have no tolerance applied to them use SPICE tolerances (see [section on SPICE tolerances](#default-tolerance-handling-spice-mode)). 

# Printing the tolerances and the natures assigned to unknowns/residuals
To print the tolerances and natures that have been assigned to individual unknowns/residuals, add the following command to the control block. 
```
print tolerances
```

If a tolerance comes from the flow nature of a `<discipline>` its name is `<discipline>.flow`. If it comes from the potential nature of a discipline, its name if `<discipline>.potential`. The nature name is not resolved all the way to the actual nature name because disciplines can override nature attributes (and along with that the value of `abstol`). The natures of integral quantities (idt natures) are resolved to their actual nature names because Verilog-A does not allow the idt nature of a potential or a flow to be overridden in discipline declaration. 
