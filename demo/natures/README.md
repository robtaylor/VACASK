# Tolerance handling in VACASK

# Nature and discipline import from .osdi files
Each .osdi file has its own set of natures and disciplines. They apply only to the moduels defined in that particular file. When computing global maxima (across the same nature) in the nonlinear solver the name of the nature is the one that is used to identify the same natures across .osdi files. You can list the natures and disciplined along with all their attributes by adding 
```
print device_file("<file name substring>")
```
This will print the details of all files whose canonical name contains the specified substring. 

```
print device_files
```
lists the canonical names of all loaded .osdi files. 

# Default tolerance handling (SPICE mode)
This is the default. You can explicitly enable it by adding 
```
options tolmode="spice"
```
to your control block. The following absolute tolerances are applied in this mode. 

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

Implicit equations (the corresponding unknowns and residuals) are treated as potential nodes. Natures whose name starts with a dot are builtin SPICE natures. Idt quantities are applied to reactive residual contributions. Currently they are used only in the element evaluation bypass algorithm. 

# Verilog-A tolerance handling
To enable it add 
```
options tolmode="va"
```
or
```
options tolmode="mixed"
```
to your control block. The `va` mode considers only absolute tolerances defined in .osdi files. Each node of every device applies an absolute tolerance to the circuit node to which it is connected. If two devices want to apply different absolute tolerances to the same node the lower value is used along with the nature from where this tolerance was taken. VACASK performs no nature compatibility checks. 

Tolerances that are applied to circuit nodes by the nodes of a device can be listed if the particular device is printed with 
```
print device("<device name>")
```
In the output there is a section named "Absolute tolerances in Verilog-A mode". In this section the absolute tolerances of the device's nodes are listed (the unknown, its integral, the residual, and the residual's integral). The names of the natures where the unknown and the residual tolerances are taken from are listed in the "Nodes" section. The tolerances for the integral of the unknown and the residual are taken from the idt attribute of the unknown/residual nature. If the idt attribute is missing the same tolerance is used as for the unknown/residual. 

Some nodes do not apply tolerances because OpenVAF currently does not expose their natures. These nodes are
- implicit equations created from ddt() and idt() operators
- implicit equations created manually in the module
- implicit equations of switch branches

Due to this some circuit nodes may not have any tolerance assigned. For such nodes no unknown/residual tolerance checks are performed. 

Builtin devices (independent and controlled linear voltage/current sources) apply the standard SPICE tolerances (see .voltage/.current/.flux/.charge above). 

When `tolmode="mixed"` is set a mixture of Verilog-A and SPICE tolerances is used. Circuit nodes that have no tolerance applied to them use SPICE tolerances (see [section on SPICE tolerances](#default-tolerance-handling-spice-mode)). 

# Printing the tolerances and the natures assigned to circuit nodes
To print the tolerances and natures that have been assigned to individual circuit nodes, add the following to the control block. 
```
print tolerances
```
This will print the assigned natures and tolerances for the unknown and the residual of each node. 

If node tolerance comes from the `electrical` flow nature of a discipline its name is `electrical.flow`. If it comes from the potential nature, its name if `electrical.potential`. The nature name is not resolved all the way to the actual nature name because disciplines can override nature attributes (and along with that the value of abstol). The natures of integral quantities (idt natures) are resolved to their actual nature names because Verilog-A prohibits the use of `<discipline>.flow` and `<discipline>.potential` as the name of the idt nature. 