v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 600 -400 1160 -100 {flags=graph
y1=-10e-6
y2=300e-6
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
divx=5
subdivx=1
dataset=-1
unitx=1
logx=0
logy=0
autoload=0
sweep=vds
hilight_wave=0
color=4
node=Vd.flow(br)
mode=Line
x2=1.48
x0=0
y0=-10e-6
rawfile=$netlist_dir/dc2.raw}
T {Ctrl-Click to execute launcher} 1130 30 0 0 0.3 0.3 {layer=11}
T {.save file can be created with IHP->"Create FET and BIP .save file"} 1130 150 0 0 0.3 0.3 {layer=11}
N 250 -160 250 -140 {
lab=GND}
N 250 -250 250 -220 {
lab=G}
N 380 -220 380 -160 {
lab=GND}
N 510 -220 510 -160 {
lab=GND}
N 380 -340 380 -280 {
lab=#net1}
N 510 -340 510 -280 {
lab=D}
N 380 -250 450 -250 {
lab=GND}
N 450 -250 450 -160 {
lab=GND}
N 380 -340 410 -340 {
lab=#net1}
N 470 -340 510 -340 {
lab=D}
N 250 -250 340 -250 {
lab=G}
C {devices/gnd.sym} 380 -160 0 0 {name=l1 lab=GND}
C {devices/gnd.sym} 250 -140 0 0 {name=l2 lab=GND}
C {devices/vsource.sym} 250 -190 0 0 {name=Vgs value="dc=1.2"}
C {devices/vsource.sym} 510 -250 0 0 {name=Vds value="dc=1.5"}
C {devices/gnd.sym} 510 -160 0 0 {name=l3 lab=GND}
C {devices/gnd.sym} 450 -160 0 0 {name=l4 lab=GND}
C {devices/title.sym} 160 -30 0 0 {name=l5 author="Copyright 2023 IHP PDK Authors"}
C {devices/ammeter.sym} 440 -340 1 0 {name=Vd}
C {lab_pin.sym} 250 -250 0 0 {name=p1 sig_type=std_logic lab=G}
C {lab_pin.sym} 510 -340 0 1 {name=p2 sig_type=std_logic lab=D}
C {sg13g2_pr/sg13_lv_nmos.sym} 360 -250 0 0 {name=M1
l=0.45u
w=1.0u
ng=1
m=1
model=sg13_lv_nmos
spiceprefix=X
}
C {sg13g2_pr/annotate_fet_params.sym} 110 -400 0 0 {name=annot1 ref=M1}
C {devices/launcher.sym} 1190 100 0 0 {name=h1
descr="OP annotate" 
tclcommand="xschem annotate_op"
}
C {devices/launcher.sym} 1190 130 0 0 {name=h2
descr="Load waves" 
tclcommand="
xschem raw_read $netlist_dir/dc2.raw
xschem setprop rect 2 0 fullxzoom
xschem setprop rect 2 0 x1 0
xschem setprop rect 2 0 x2 1.48
xschem setprop rect 2 0 y1 -10e-6
xschem setprop rect 2 0 y2 300e-6
"
}
C {launcher.sym} 1190 70 0 0 {name=h3
descr=SimulateVACASK
tclcommand="
# Setup the default simulation commands if not already set up
# for example by already launched simulations.
set_sim_defaults
puts $sim(spectre,0,cmd) 

# change the simulator to be used (#0 in spectre category is VACASK)
set sim(spectre,default) 0

# Create FET and BIP .save file
mkdir -p $netlist_dir
write_data [save_params] $netlist_dir/[file rootname [file tail [xschem get current_name]]].save

# run netlist and simulation
xschem netlist
simulate
"}
C {simulator_commands_shown.sym} 60 80 0 0 {
name=Libs_VACASK
simulator=vacask
only_toplevel=false
value="
include \\"sg13g2_vacask_common.lib\\"
include \\"cornerMOSlv.lib\\" section=mos_tt
include \\"cornerMOShv.lib\\" section=mos_tt
include \\"cornerHBT.lib\\" section=hbt_typ
include \\"cornerRES.lib\\" section=res_typ
include \\"cornerCAP.lib\\" section=cap_typ
"
      }
C {simulator_commands_shown.sym} 420 80 0 0 {
name=Script_VACASK
simulator=vacask
only_toplevel=false
value="
// Use single quotes where possible (Python)
// If double quotes are needed, escape them. 

control
  // In order to display operating point the 
  // raw file name must be the same as the netlist name

  // Save operating point data
  include \\"dc_lv_nmos.save\\"
  save default

  analysis dc_lv_nmos op

  sweep vds instance=\\"Vds\\" parameter=\\"dc\\" from=0 to=1.5 step=0.01
    analysis dc1 op
  
  sweep vgs instance=\\"Vgs\\" parameter=\\"dc\\" from=0.2 to=1.4 step=0.2
  sweep vds instance=\\"Vds\\" parameter=\\"dc\\" from=0 to=1.5 step=0.01
    analysis dc2 op
  
  postprocess(PYTHON, \\"runme.py\\")
endc

embed \\"runme.py\\" <<<FILE
from rawfile import rawread
import numpy as np
import matplotlib.pyplot as plt

fig1 = plt.figure(figsize=(6,4), dpi=100)
fig1.suptitle('LV NMOS')
ax_dict = fig1.subplot_mosaic('A')

ax_dict['A'].set_ylabel('Id [uA]')
ax_dict['A'].set_xlabel('Vds [V]')

raw = rawread('dc2.raw')
plot = raw.get(sweeps=1)
for ii in range(plot.sweepGroups):
    sdata = plot.sweepData(ii)
    traceName = 'vgs=%.2e' % sdata['vgs']
    ax_dict['A'].plot(
        plot[ii, 'D'], -plot[ii, 'Vds:flow(br)']*1e6, 
        label=traceName
    )
ax_dict['A'].legend(loc='upper right')
fig1.tight_layout()
plt.show()
>>>FILE
"
      }
