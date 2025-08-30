v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 -160 -630 640 -230 {flags=graph
y1=-0.1
y2=1
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=10e-9
divx=5
subdivx=1


dataset=-1
unitx=1
logx=0
logy=0
color="4 5"
node="out
in"}
N -180 140 -180 160 {
lab=GND}
N -180 20 -180 80 {
lab=in}
N -180 20 -140 20 {
lab=in}
N -200 20 -180 20 {
lab=in}
N 0 20 60 20 {
lab=out}
N 0 20 0 50 {lab=out}
N -60 20 0 20 {
lab=out}
N -180 -190 -180 -160 {lab=VDD}
N -180 -100 -180 -80 {lab=GND}
N -180 -190 -140 -190 {lab=VDD}
N -50 -190 -50 -160 {lab=VSS}
N -50 -100 -50 -80 {lab=GND}
N -50 -190 -10 -190 {lab=VSS}
C {devices/gnd.sym} -180 160 0 0 {name=l2 lab=GND}
C {devices/vsource.sym} -180 110 0 0 {name=Vin value="type=\\"pulse\\" val0=0 val1=0.9 delay=0 rise=100p fall=100p width=2n period=4n"}
C {devices/title.sym} -130 260 0 0 {name=l5 author="Copyright 2024 IHP PDK Authors"}
C {devices/launcher.sym} 230 -170 0 0 {name=h5
descr="load waves Ctrl + left click" 
tclcommand="xschem raw_read $netlist_dir/tran1.raw
xschem setprop rect 2 0 x1 0
xschem setprop rect 2 0 x2 10e-9
xschem setprop rect 2 0 y1 -0.1
xschem setprop rect 2 0 y2 1
"
}
C {devices/lab_pin.sym} -200 20 0 0 {name=p1 sig_type=std_logic lab=in}
C {devices/lab_pin.sym} 60 20 2 0 {name=p2 sig_type=std_logic lab=out}
C {sg13g2_stdcells/sg13g2_inv_1.sym} -100 20 0 0 {name=x1 VDD=VDD VSS=VSS prefix=sg13g2_ }
C {devices/parax_cap.sym} 0 60 0 0 {name=C2 gnd=GND value=4f m=1}
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
C {simulator_commands_shown.sym} 560 -90 0 0 {
name=VACASK_SIMULATION
simulator=vacask
only_toplevel=false
value="
control
  options temp=27
  analysis tran1 tran step=50p stop=20n
endc
"
      }
C {devices/vsource.sym} -180 -130 0 0 {name=Vdd value="dc=0.9"}
C {devices/gnd.sym} -180 -80 0 0 {name=l1 lab=GND}
C {devices/lab_pin.sym} -140 -190 2 0 {name=p3 sig_type=std_logic lab=VDD}
C {devices/vsource.sym} -50 -130 0 0 {name=Vss value="dc=0"}
C {devices/gnd.sym} -50 -80 0 0 {name=Vss1 lab=GND
value="dc=0"}
C {devices/lab_pin.sym} -10 -190 2 0 {name=p4 sig_type=std_logic lab=VSS}
