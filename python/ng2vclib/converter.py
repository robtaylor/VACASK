from .m_file import FileLoaderMixin
from .m_masters import MastersMixin
from .m_output import OutputMixin
from .m_model import ModelMixin
from .m_inst_n import InstanceNMixin

class Converter(
    FileLoaderMixin, 
    MastersMixin, 
    OutputMixin, 
    ModelMixin, 
    InstanceNMixin
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
    * as_toplevel .. treat input file as toiplevel netlist. Possible values 
      * auto .. detect automatically based on .end
      * no .. treat as toplevel
      * yes .. do not treat as toplevel
    * columns .. number of columns per line
    * flat .. produce flat netlist (resolve all .incliues and .libs)

    The data member dictionary holds the following keys:
    * title .. unprocessed first line of toplevel file
    * control .. control block lines
    * models .. a dictionary of models, key is subcircuit name
      where None represents the toplevel circuit. values are dictionaries
      with model name as key holding tuples of the form
      holding elements:
      * builtin .. is it a SPICE builtin
      * model type (npn, pnp, nmos, pmos, ...)
      * level (integer)
      * version (string)
      * params .. list of (name, value) pairs. Parameters level and version are 
        removed from this list if they had to be collected. 
    * subckts .. a dictionary of subcircuit definitions with name for key
      and values that are tuples holding
      * terminals .. list of terminals
      * parameters .. list of (name, value) pairs for parameters
    * is_toplevel .. True if input file is a toplevel netlist
    """
    def __init__(self, cfg):
        self.cfg = cfg
        self.data = {
            "control": [], 
            "models": {}, 
            "subckts": {}, 
        }
    
    def read(self, filename):
        deck = self.read_file(filename)
        self.data["deck"] = deck
        return deck

