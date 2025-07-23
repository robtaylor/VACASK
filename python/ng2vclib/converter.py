from .m_file import FileLoaderMixin
from .m_masters import MastersMixin
from .m_output import OutputMixin
from .m_model import ModelMixin
from .m_inst_passive import InstancePassiveMixin
from .m_inst_d import InstanceDMixin
from .m_inst_n import InstanceNMixin
from .m_inst_q import InstanceQMixin
from .m_devices import DevicesMixin

from . import dfl

class Converter(
    FileLoaderMixin, 
    MastersMixin, 
    OutputMixin, 
    ModelMixin, 
    InstancePassiveMixin, 
    InstanceDMixin, 
    InstanceNMixin, 
    InstanceQMixin, 
    DevicesMixin
):
    """
    Ngspice to VACASK netlist converter. 
    *cfg* is the configuration dictionary with the following keys: 
    * sourcepath .. a list of directories where .include and .lib files
      are located
    * type_map .. a list of model types (tuples) with entries
      * SPICE name (npn, pnp, nmos, pmos, ...)
      * parameters to add when converting to VACASK
      * family (used as key in other tables)
      * remove_level .. True if level parameter should be removed 
      * remove_version .. True if version parameter should be removed
    * family_map .. maps tuples of the form
      * family
      * level (integer, None is default)
      * version (string, None is default)
      into tuples of the form
      * osdi file to load
      * module name
    * default_models .. maps device letters to default models
      Values are tuples holding
      * osdi file name
      * module name
    * default_model_prefix .. prefix for default model names
    * as_toplevel .. treat input file as toiplevel netlist. Possible values 
      * auto .. detect automatically based on .end
      * no .. treat as toplevel
      * yes .. do not treat as toplevel
    * all_models .. if True writes all models to the output, 
      otherwise ony used models are written
    * read_depth .. include/lib depth up to which the files are read.
      Default (None) is equivalent to no depth limit. 
    * process_depth .. include/lib depth up to which the files are processed. 
      Default (None) is equivalent to no depth limit. 
    * output_depth .. include/lib depth up to which the files are written
      to a VACASK netlist. Default (None) is equivalent to no depth limit. 
    * patch .. a dictionary of patches with file name for key. Files are 
      matched to this key with the last part of their name. 
      Values are lists of pairs of the form (old, new). If the old string is 
      found in the beginning of a line the whole line is replaced with the 
      new string. 
    * columns .. number of columns per line
    
    The data member dictionary holds the following keys:
    * title .. unprocessed first line of toplevel file
    * control .. control block lines
    * models .. a dictionary of available models, key is subcircuit name
      where None represents the toplevel circuit. values are dictionaries
      with model name as key holding tuples of the form
      holding elements:
      * builtin .. is it a SPICE builtin
      * model type (npn, pnp, nmos, pmos, ...)
      * model family
      * level (integer)
      * version (string)
      * params .. list of (name, value) pairs. 
    * model_usage .. a dictionary tracking where individual models are used. 
      Key is a pair of the form (model name, subcircuit where model is defined). 
      Value is a set of subcircuit names (None for toplevel circuit) 
      where this model is used. 
    * default_models_needed .. set of device letters
    * osdi_loads .. set of files loaded by pre_osdi in control block
    * subckts .. a dictionary of subcircuit definitions with name for key
      and values that are tuples holding
      * terminals .. list of terminals
      * parameters .. list of (name, value) pairs for parameters
    * is_toplevel .. True if input file is a toplevel netlist
    """
    def __init__(self, cfg=None):
        # If no config is given, use default config
        if cfg is None:
          cfg = dfl.default_config()
        self.cfg = cfg
        self.data = {
            "control": [], 
            "models": {}, 
            "subckts": {}, 
            "model_usage": {}, 
            "default_models_needed": set(), 
            "osdi_loads": set(), 
        }
    
    def read(self, filename):
      """
      Read a file, store deck. 
      """
      deck = self.read_file(filename)
      self.data["deck"] = deck
      return deck

    def convert(self, fromFile, toFile=None):
      self.read(fromFile)
      
      self.collect_masters()

      out = self.vacask_file()

      if toFile is None:
        for l in out:
          print(l)
      else:
        with open(toFile, "w") as f:
          for l in out:
            f.write(l)
            f.write("\n")

