// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { INTERPROSCAN_MERGE } from '/dados04/matheus/rgaprofiler/omatheuspimenta-rgaprofiler/modules/local/interproscan_merge/tests/../main.nf'

workflow {

    // define custom rules for JSON that will be generated.
    def jsonOutput = createJsonOutput()
    def jsonWorkflowOutput = createJsonWorkflowOutput()

    def input = []

    // run dependencies
    

    // process mapping
    input = []
    
                def chunk1 = file("${workDir}/chunk1_interpro.tsv")
                chunk1.text = "SoffiXsponR570.7os1g018900.1.p\tMD5A\t76\tPfam\tPF00001\tDomain one\t1\t28\t1.0E-5\tT\t01-01-2026\t-\t-\t-\t-\n"
                def chunk2 = file("${workDir}/chunk2_interpro.tsv")
                chunk2.text = "SoffiXsponR570.10Eg022100.1.p\tMD5B\t1160\tPfam\tPF00931\tNB-ARC domain\t202\t371\t3.1E-6\tT\t01-01-2026\tIPR002182\tNB-ARC\t-\t-\n"

                input[0] = [
                    [ id:'r570_subsample', single_end:false ], // meta
                    [ chunk1, chunk2 ]
                ]
                
    //----

    //run process
    INTERPROSCAN_MERGE.run(input.toArray())

    if (INTERPROSCAN_MERGE.output){

        // consumes all named output channels and stores items in a json file
        INTERPROSCAN_MERGE.out.getNames().each { name ->
            serializeChannel(name, INTERPROSCAN_MERGE.out.getProperty(name), jsonOutput, params.nf_test_output)
        }	  

        // consumes all unnamed output channels and stores items in a json file
        def array = INTERPROSCAN_MERGE.out as List<Object>
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