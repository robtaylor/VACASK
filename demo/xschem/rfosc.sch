v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 450 -180 1250 220 {flags=graph
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
divx=5
subdivx=1
xlabmag=1.0
ylabmag=1.0
node=c
color=4
dataset=-1
unitx=1
logx=0
logy=0
rawfile=$netlist_dir/tran1.raw
x2=10e-6
x1=0
y1=-2
y2=25}
N 40 -100 40 -50 {lab=vcc}
N 40 -100 370 -100 {lab=vcc}
N 370 -100 370 20 {lab=vcc}
N 150 -100 150 -80 {lab=vcc}
N 230 -100 230 -80 {lab=vcc}
N 150 -20 150 20 {lab=c}
N 230 -20 230 20 {lab=c}
N 150 0 230 0 {lab=c}
N 150 80 150 130 {lab=e}
N 150 110 230 110 {lab=e}
N 230 80 230 110 {lab=e}
N 40 10 40 100 {lab=b}
N 40 50 110 50 {lab=b}
N 150 190 150 210 {lab=GND}
N 40 210 150 210 {lab=GND}
N 40 160 40 210 {lab=GND}
N 150 210 370 210 {lab=GND}
N 370 80 370 210 {lab=GND}
N 150 210 150 220 {lab=GND}
C {npn.sym} 130 50 0 0 {name=Q1
model=q2n2222
device=2n2222
footprint=SOT23
area=1
m=1}
C {res.sym} 40 -20 0 0 {name=R1
value=47k
footprint=1206
device=resistor
m=1}
C {res.sym} 150 160 0 0 {name=R2
value=470
footprint=1206
device=resistor
m=1}
C {capa.sym} 150 -50 0 0 {name=C1
m=1
value=20p
footprint=1206
device="ceramic capacitor"}
C {capa.sym} 230 50 0 0 {name=C2
m=1
value=10p
footprint=1206
device="ceramic capacitor"}
C {capa.sym} 40 130 0 0 {name=C3
m=1
value=1n
footprint=1206
device="ceramic capacitor"}
C {ind.sym} 230 -50 0 0 {name=L1
m=1
value=100n
footprint=1206
device=inductor}
C {vsource.sym} 370 50 0 0 {name=Vcc value="dc=9" savecurrent=false}
C {gnd.sym} 150 220 0 0 {name=l2 lab=GND}
C {lab_wire.sym} 40 50 0 0 {name=p1 sig_type=std_logic lab=b
}
C {lab_wire.sym} 150 0 0 0 {name=p2 sig_type=std_logic lab=c
}
C {lab_wire.sym} 150 110 0 0 {name=p3 sig_type=std_logic lab=e
}
C {lab_wire.sym} 190 -100 0 0 {name=p4 sig_type=std_logic lab=vcc
}
C {simulator_commands_shown.sym} 20 330 0 0 {name=Models
simulator=vacask
only_toplevel=false 
value="
load \\"spice/bjt.osdi\\"
model q2n2222 sp_bjt ( type=1
  is=15.2f nf=1 bf=105 vaf=98.5 ikf=0.5
  ise=8.2p ne=2 br=4 nr=1 var=20 ikr=0.225
  re=0.373 rb=1.49 rc=0.149 xtb=1.5
  cje=35.5p cjc=12.2p tf=500p tr=85n
)
"}
C {simulator_commands_shown.sym} 500 330 0 0 {name=Commands
simulator=vacask
only_toplevel=false 
value="
control
  options reltol=1e-6
  analysis op1 op
  analysis tran1 tran step=0.1n stop=10u maxstep=0.1n icmode=\\"uic\\"
endc
"}
