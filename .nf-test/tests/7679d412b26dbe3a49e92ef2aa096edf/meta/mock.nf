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
    
                // Extract just the RGA sequence out of the shared subsample FASTA (plain
                // Groovy string parsing -- splitFasta is a channel operator, not usable
                // as plain script code in this synchronous test block).
                def subsample_records = file("${projectDir}/tests/testdata/r570_subsample.fasta").text.split('(?=\n>)')
                def rga_record = subsample_records.find { it.contains('SoffiXsponR570.10Eg022100.1.p') }
                def rga_fasta = file("${workDir}/rga.fasta")
                rga_fasta.text = rga_record.trim() + '\n'

                input[0] = [
                    [ id:'rga_test', single_end:false ],
                    rga_fasta
                ]
                input[1] = file("${projectDir}/softwares/InterProScan/interproscan-5.78-109.0", checkIfExists: true)
                
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