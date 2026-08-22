// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { INTERPROSCAN } from '/dados04/matheus/rgaprofiler/omatheuspimenta-rgaprofiler/modules/local/interproscan/tests/../main.nf'

workflow {

    // define custom rules for JSON that will be generated.
    def jsonOutput = createJsonOutput()
    def jsonWorkflowOutput = createJsonWorkflowOutput()

    def input = []

    // run dependencies
    

    // process mapping
    input = []
    
                // Create fasta
                def chunk_fasta = file("${workDir}/chunk_1.fasta")
                chunk_fasta.text = '>Seq1\nMKVLLLLSAVIGASACQGGVAAQNCYKQRCARRCGTRCGLCCQAQCAKACGCCCA\n'

                //Create a "Fake" folder of InterProScan
                def mock_ipr_dir = file("${workDir}/fake_interproscan_5.78")
                mock_ipr_dir.mkdirs()
                
                def mock_sh = file("${mock_ipr_dir}/interproscan.sh")
                mock_sh.text = '''#!/bin/bash
                if [[ "$1" == "--version" ]]; then
                    echo "InterProScan version 5.78-109.0"
                    exit 0
                fi
                # Simulate the expected output
                echo -e "Seq1\t...fake_interpro_results..." > teste_ipr_interpro.tsv
                '''.stripIndent()
                
                // mock permission
                mock_sh.setPermissions('rwxr-xr-x')

                // fed the chanels
                input[0] = [
                    [ id:'teste_ipr', single_end:false ],
                    chunk_fasta
                ]
                input[1] = mock_ipr_dir // mock folder as argument
                
    //----

    //run process
    INTERPROSCAN.run(input.toArray())

    if (INTERPROSCAN.output){

        // consumes all named output channels and stores items in a json file
        INTERPROSCAN.out.getNames().each { name ->
            serializeChannel(name, INTERPROSCAN.out.getProperty(name), jsonOutput, params.nf_test_output)
        }	  

        // consumes all unnamed output channels and stores items in a json file
        def array = INTERPROSCAN.out as List<Object>
        def i = 0
        array.each { output ->
            serializeChannel(i, output, jsonOutput, params.nf_test_output)
            i += 1
        }    	

    }

    // get topics

    // finalize test
    workflow.onComplete = {
        def result = [
            success: workflow.success,
            exitStatus: workflow.exitStatus,
            errorMessage: workflow.errorMessage,
            errorReport: workflow.errorReport
        ]
        new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
        
    }
}

def serializeChannel(name, channel, jsonOutput, outputDir) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: { entry ->
            list.add(entry)
        },
        onComplete: {
            def map = new HashMap()
            map[_name] = list
            def filename = "${outputDir}/output_${_name}.json"
            new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}

def serializeTopic(name, topic, jsonOutput, outputDir) {
    def list = [ ]
    topic.subscribe(
        onNext: { entry ->
            list.add(entry)
        },
        onComplete: {
            def map = new HashMap()
            map[name] = list
            def filename = "${outputDir}/topic_${name}.json"
            new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}

def createJsonOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the process
    return [
        toJson: { obj ->
            def converted = convertPathsToStrings(obj)
            return groovy.json.JsonOutput.toJson(converted)
        }
    ]
}

def convertPathsToStrings(obj) {
    if (obj instanceof java.nio.file.Path) {
        return obj.toAbsolutePath().toString()
    } else if (obj instanceof Map) {
        return obj.collectEntries { k, v -> [k, convertPathsToStrings(v)] }
    } else if (obj instanceof Collection) {
        return obj.collect { it -> convertPathsToStrings(it) }
    } else {
        return obj
    }
}

def createJsonWorkflowOutput(_input = null) {
    // _input is needed because a closure is provided to all functions called in the workflow
    return [
        toJson: { obj ->
            def filtered = removeNullValues(obj)
            return groovy.json.JsonOutput.toJson(filtered)
        }
    ]
}

def removeNullValues(obj) {
    if (obj instanceof Map) {
        return obj.findAll { _k, v -> v != null }.collectEntries { k, v -> [k, removeNullValues(v)] }
    } else if (obj instanceof Collection) {
        return obj.findAll { it -> it != null }.collect { it -> removeNullValues(it) }
    } else {
        return obj
    }
}