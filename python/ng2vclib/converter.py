from .m_file import FileLoaderMixin
from .m_masters import MastersMixin

class Converter(
    FileLoaderMixin, 
    MastersMixin
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
      * has_level .. True if model has level parameter that needs to be collected
      * has_version .. True if model has version parameter that needs to be collected
    * family_map .. maps tuples of the form
      * family
      * level (integer, None is default)
      * version (string, None is default)
      into tuples of the form
      * osdi file to load
      * module name
    
    The data member dictionary holds the following keys:
    * lines .. file lines without control block
    * control .. control block lines
    * models .. a dictionary of models, key is subcircuit name
      where None represents the toplevel circuit. values are dictionaries
      with model name as key holding tuples of the form
      holding elements:
      * depth .. inclusion depth (toplevel in 0)
      * builtin .. is it a SPICE builtin
      * model type (npn, pnp, nmos, pmos, ...)
      * level (integer)
      * version (string)
      * params .. list of (name, value) pairs. Parameters level and version are 
        removed from this list if they had to be collected. 
    * subckts .. a dictionary of subcircuit definitions with name for key
      and values that are tuples holding
      * depth .. inclusion depth (toplevel in 0)
      * terminals
      * parameters
    """
    def __init__(self, cfg):
        self.cfg = cfg
        self.data = {
            "control": [], 
            "models": {}, 
            "subckts": {}, 
        }
    
