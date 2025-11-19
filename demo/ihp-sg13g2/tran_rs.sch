v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 -160 -630 640 -230 {flags=graph
y1=-0.1
y2=8
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=14e-9
divx=5
subdivx=1


dataset=-1
unitx=1
logx=0
logy=0
color="4 4 10 10"
node="\\"s;s 6 +\\"
\\"r;r 4 +\\"
\\"q; q 2 +\\"
qbar"
digital=0}
N -180 -190 -180 -160 {lab=VDD}
N -180 -100 -180 -80 {lab=GND}
N -180 -190 -140 -190 {lab=VDD}
N -50 -190 -50 -160 {lab=VSS}
N -50 -100 -50 -80 {lab=GND}
N -50 -190 -10 -190 {lab=VSS}
N 270 180 390 210 {lab=q}
N 270 210 390 180 {lab=qbar}
N 270 140 270 180 {lab=q}
N 270 210 270 250 {lab=qbar}
N 390 120 390 180 {lab=qbar}
N 390 210 390 270 {lab=q}
N 390 120 440 120 {lab=qbar}
N 390 270 440 270 {lab=q}
N 170 100 270 100 {lab=s}
N 170 290 270 290 {lab=r}
C {devices/launcher.sym} 230 -170 0 0 {name=h5
descr="load waves Ctrl + left click" 
tclcommand="xschem raw_read $netlist_dir/tran1.raw
xschem setprop rect 2 0 x1 0
xschem setprop rect 2 0 x2 14e-9
xschem setprop rect 2 0 y1 -0.1
xschem setprop rect 2 0 y2 8
"
}
C {simulator_commands_shown.sym} 150 -90 0 0 {
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
include \\"sg13g2_stdcell.inc\\"
"
      }
C {devices/vsource.sym} -180 -130 0 0 {name=Vdd value="dc=0.9"}
C {devices/gnd.sym} -180 -80 0 0 {name=l1 lab=GND}
C {devices/lab_pin.sym} -140 -190 2 0 {name=p3 sig_type=std_logic lab=VDD}
C {devices/vsource.sym} -50 -130 0 0 {name=Vss value="dc=0"}
C {devices/gnd.sym} -50 -80 0 0 {name=Vss1 lab=GND
value="dc=0"}
C {devices/lab_pin.sym} -10 -190 2 0 {name=p4 sig_type=std_logic lab=VSS}
C {devices/lab_pin.sym} 170 100 0 0 {name=p5 sig_type=std_logic lab=s}
C {devices/lab_pin.sym} 170 290 0 0 {name=p6 sig_type=std_logic lab=r}
C {devices/lab_pin.sym} 440 120 0 1 {name=p9 sig_type=std_logic lab=qbar}
C {devices/lab_pin.sym} 440 270 0 1 {name=p10 sig_type=std_logic lab=q}
C {devices/vsource.sym} 220 130 0 1 {name=Vin1 value="type=\\"pulse\\" val0=0 val1=0.9 delay=1n rise=100p fall=100p width=2n period=9n"}
C {devices/vsource.sym} 220 320 0 1 {name=Vin2 value="type=\\"pulse\\" val0=0 val1=0.9 delay=5n rise=100p fall=100p width=2n"}
C {devices/gnd.sym} 220 160 0 0 {name=l2 lab=GND}
C {devices/gnd.sym} 220 350 0 0 {name=l3 lab=GND}
C {sg13g2_stdcells/sg13g2_nor2_1.sym} 330 120 0 0 {name=x1 VDD=VDD VSS=VSS prefix=sg13g2_ }
C {sg13g2_stdcells/sg13g2_nor2_1.sym} 330 270 0 0 {name=x2 VDD=VDD VSS=VSS prefix=sg13g2_ }
C {command_block.sym} 590 -140 0 0 {name=CMD
only_toplevel=false
}
C {tran.sym} 590 120 0 0 {name=tran1
only_toplevel=false 
order="1"
sweep=""
step="50p"
stop="14n"
start=""
maxstep=""
icmode=""
nodeset=""
ic=""
store=""
write=""
}
C {verbatim.sym} 590 -20 0 0 {name=verbatim1
only_toplevel=false 
order="0"
simulator="vacask"
verbatim="options temp=27"
}
