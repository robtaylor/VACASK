set SG13G2_MODELS_VACASK ${PDK_ROOT}/${PDK}/libs.tech/vacask/models
puts stderr "SG13G2_MODELS_VACASK: $SG13G2_MODELS_VACASK"

# Rename write_save_lines for spice
rename write_save_lines write_save_lines_spice

# Define write_save_lines for vacask
proc write_save_lines_vacask {type model schpath spiceprefix instname} {
  global sch_expand
  if {[regexp {[pn]mos} $type]} {
    set m n$model
    set devpath [ string map {. :} '$schpath$spiceprefix$instname.$m' ]

    append sch_expand(savelist) "save p($devpath,ids)\n"
    append sch_expand(savelist) "save p($devpath,gm)\n"
    append sch_expand(savelist) "save p($devpath,gds)\n"
    append sch_expand(savelist) "save p($devpath,vth)\n"
    append sch_expand(savelist) "save p($devpath,vgs)\n"
    append sch_expand(savelist) "save p($devpath,vdss)\n"
    append sch_expand(savelist) "save p($devpath,vds)\n"
    append sch_expand(savelist) "save p($devpath,cgg)\n"
    append sch_expand(savelist) "save p($devpath,cgsol)\n"
    append sch_expand(savelist) "save p($devpath,cgdol)\n"
  } elseif {[regexp {vertical_npn} $type]} {
    if {[regexp {_5t$} $model]} {
      set model [string range $model 0 end-3]
    }
    set m q$model
    set devpath [ string map {. :} '$schpath$spiceprefix$instname.$m' ]

    # TODO: Verilog-A VBIC does not expose output variables
    # append sch_expand(savelist) "save p($devpath,gm)\n"
    # append sch_expand(savelist) "save p($devpath,go)\n"
    # append sch_expand(savelist) "save p($devpath,gmu)\n"
    # append sch_expand(savelist) "save p($devpath,gpi)\n"
    # append sch_expand(savelist) "save p($devpath,gx)\n"
    # append sch_expand(savelist) "save p($devpath,vbe)\n"
    # append sch_expand(savelist) "save p($devpath,vbc)\n"
    # append sch_expand(savelist) "save p($devpath,ib)\n"
    # append sch_expand(savelist) "save p($devpath,ic)\n"
    # append sch_expand(savelist) "save p($devpath,cbe)\n"
    # append sch_expand(savelist) "save p($devpath,cbc)\n"
    # append sch_expand(savelist) "save p($devpath,cbep)\n"
    # append sch_expand(savelist) "save p($devpath,cbcp)\n"
  }
}

# Common wrapper
proc write_save_lines {args} {
  if [ sim_is_vacask ] {
     return [ write_save_lines_vacask {*}$args ]
  } else {
     return [ write_save_lines_spice {*}$args ]
  }
}

# Rename save_params for spice
rename save_params save_params_spice

# Common wrapper
proc save_params {args} {
  if [ sim_is_vacask ] {
     return "// [ save_params_spice {*}$args ]"
  } else {
     return [ save_params_spice {*}$args ]
  }
}

# Rename MOSFET annotator
rename display_fet_params display_fet_params_spice

# MOSFET annotator for VACASK
proc display_fet_params_vacask {instname} {
  set txt {}
  set schpath [xschem get sim_sch_path]
  set symbol [xschem getprop instance $instname cell::name]
  set spiceprefix [xschem getprop instance $instname spiceprefix]
  set model [xschem translate $instname @model]
  set type [xschem getprop symbol $symbol type]

  if {[regexp {[pn]mos} $type]} {
    set m n$model
    # Do not replace . with : because raw file variables are renamed
    # when raw file is loaded so that . is used as the separator
    set devpath $schpath$spiceprefix$instname.$m

    append txt "ids   = [to_eng [xschem raw value $devpath\.ids -1]]\n"
    append txt "gm    = [to_eng [xschem raw value $devpath\.gm -1]]\n"
    append txt "gds   = [to_eng [xschem raw value $devpath\.gds -1]]\n"
    append txt "vth   = [to_eng [xschem raw value $devpath\.vth -1]]\n"
    append txt "vgs   = [to_eng [xschem raw value $devpath\.vgs -1]]\n"
    append txt "vdss  = [to_eng [xschem raw value $devpath\.vdss -1]]\n"
    append txt "vds   = [to_eng [xschem raw value $devpath\.vds -1]]\n"
    append txt "cgg   = [to_eng [xschem raw value $devpath\.cgg -1]]\n"
    set pi 3.141592654
    set gm [xschem raw value $devpath\.gm -1]
    set cgg [xschem raw value $devpath\.cgg -1]
    set cgdol [xschem raw value $devpath\.cgdol -1]
    set cgsol [xschem raw value $devpath\.cgsol -1]
    set ids [xschem raw value $devpath\.ids -1]
    if {[catch { expr $gm / $ids} gmid]} {
      set gmid {}
    }
    if {[catch { expr $gm / 2 / $pi / ($cgg + $cgdol + $cgsol)} ft]} {
      set ft {}
    }
    append txt "ft    = [to_eng ${ft}]\n"
    append txt "gm/id = [to_eng [expr $gmid]]\n"
  }
  return $txt
}

# Common wrapper
proc display_fet_params {args} {
  if [ sim_is_vacask ] {
     return [ display_fet_params_vacask {*}$args ]
  } else {
     return [ display_fet_params_spice {*}$args ]
  }
}

# TODO: annotator for HBT

# Menu manipulation
proc menuvacask {} {
  global has_x netlist_dir
  if { [info exists has_x] } {
    set topwin [xschem get top_path]

    # VACASK libraries
    $topwin.menubar.ihp insert 3 command -label {Add VACASK models symbol} -command {
      xschem place_symbol devices/simulator_commands_shown.sym {
name=Libs_VACASK
simulator=vacask
only_toplevel=false
value="
include \"sg13g2_vacask_common.lib\"
include \"cornerMOSlv.lib\" section=mos_tt
include \"cornerMOShv.lib\" section=mos_tt
include \"cornerHBT.lib\" section=hbt_typ
include \"cornerRES.lib\" section=res_typ
include \"cornerCAP.lib\" section=cap_typ
"
      }
    }
  }
}

append postinit_commands "menuvacask\n"
