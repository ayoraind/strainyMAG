include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utilities.nf'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:] as nextflow.script.ScriptBinding$ParamsMap
    // Defaults for configurable variables
    params.help = false
    params.version = false
    params.reads = false
    params.valid_mode = false
    params.output_dir = false
    return params
}

def check_params(Map params) { 
    final_params = params
    
    // set up reads files
    final_params.reads = check_mandatory_parameter(params, 'reads')
     
    // set up output directory
    final_params.output_dir = check_mandatory_parameter(params, 'output_dir') - ~/\/$/
    
    // check if valid mode is valid
    final_params.valid_mode = check_parameter_value('mode', final_params.valid_mode, ['--pacbio-raw', '--pacbio-corr', '--pacbio-hifi', '--nano-raw', '--nano-corr', '--nano-hq'])
         	
    return final_params
}

