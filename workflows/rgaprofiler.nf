/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { logColours             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rgaprofiler_pipeline'

include { FASTA_QC           } from '../modules/local/fasta_qc'
include { DEEPCOIL2          } from '../modules/local/deepcoil2'
include { PHOBIUS            } from '../modules/local/phobius'
include { INTERPROSCAN       } from '../modules/local/interproscan'
include { INTERPROSCAN_MERGE } from '../modules/local/interproscan_merge'
include { DEEPLOC2           } from '../modules/local/deeploc2'
include { SIGNALP6           } from '../modules/local/signalp6'
include { DEEPTMHMM          } from '../modules/local/deeptmhmm'
include { RGA_CLASSIFY       } from '../modules/local/rga_classify'
include { RGA_REPORT         } from '../modules/local/rga_report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RGAPROFILER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir

    main:

    // Check the InterProScan database dir path exists and has the files this pipeline expects
    // (see bin/check_software_present.sh and docs/software-setup.md).
    if (!params.interproscan_db) {
        exit 1, "ERROR: Path to InterProScan database not provided! Use --interproscan_db /path/to/directory"
    }
    def ipr_dir = file(params.interproscan_db, checkIfExists: true)
    def ipr_check = ["${projectDir}/bin/check_software_present.sh", "interproscan", ipr_dir.toString()].execute()
    ipr_check.waitFor()
    if (ipr_check.exitValue() != 0) {
        log.error(ipr_check.err.text)
        exit 1, "ERROR: InterProScan installation at ${ipr_dir} is incomplete. See message above."
    }

    // Check DeepLoc2's DTU model checkpoints + ESM1b base encoder cache are present
    // (see bin/check_software_present.sh and docs/software-setup.md).
    def deeploc2_check = ["${projectDir}/bin/check_software_present.sh", "deeploc2", params.softwares_dir].execute()
    deeploc2_check.waitFor()
    if (deeploc2_check.exitValue() != 0) {
        log.error(deeploc2_check.err.text)
        exit 1, "ERROR: DeepLoc2 model weights under ${params.softwares_dir}/DeepLoc2 are incomplete. See message above."
    }
    def deeploc2_models = file("${params.softwares_dir}/DeepLoc2/DeepLoc2/models", checkIfExists: true)
    def deeploc2_torch_cache = file("${params.softwares_dir}/DeepLoc2/torch_cache", checkIfExists: true)

    // SignalP6 has no runtime --device flag: GPU vs CPU is baked into the weight files
    // themselves (see modules/local/signalp6 and docs/software-setup.md), so — unlike
    // every other GPU-capable tool here, where the same staged weights work either way —
    // which *directory* of weights to stage has to be decided here, before the
    // channel/task is even built. This mirrors the exact GPU-resolution logic in
    // conf/base.config's process_gpu label (see that file's comment for why it must be
    // resolved lazily and not via a params.* value); the two must stay in sync.
    def signalp6_use_gpu
    if (params.use_gpu == 'auto') {
        def gpu_probe = ["bash", "${projectDir}/bin/detect_gpu.sh"].execute()
        signalp6_use_gpu = (gpu_probe.waitFor() == 0)
    } else {
        signalp6_use_gpu = (params.use_gpu == 'true')
    }
    def signalp6_models_subdir = signalp6_use_gpu ? 'models_gpu' : 'models'

    // Check SignalP6's DTU model weights are present (see bin/check_software_present.sh
    // and docs/software-setup.md). Passing signalp6_use_gpu makes the GPU-converted copy
    // a hard requirement (not just an optional extra) whenever GPU execution was actually
    // requested/detected, with a message telling the user exactly how to produce it.
    def signalp6_check = ["${projectDir}/bin/check_software_present.sh", "signalp6", params.softwares_dir, signalp6_use_gpu.toString()].execute()
    signalp6_check.waitFor()
    if (signalp6_check.exitValue() != 0) {
        log.error(signalp6_check.err.text)
        exit 1, "ERROR: SignalP6 model weights under ${params.softwares_dir}/SignalP6 are incomplete. See message above."
    }
    def signalp6_models = file("${params.softwares_dir}/SignalP6/signalp-6-package/${signalp6_models_subdir}", checkIfExists: true)

    // Check DeepTMHMM's DTU/BioLib model weights are present (see bin/check_software_present.sh
    // and docs/software-setup.md). CC BY-NC-SA 4.0 -- academic/non-commercial use only.
    def deeptmhmm_check = ["${projectDir}/bin/check_software_present.sh", "deeptmhmm", params.softwares_dir].execute()
    deeptmhmm_check.waitFor()
    if (deeptmhmm_check.exitValue() != 0) {
        log.error(deeptmhmm_check.err.text)
        exit 1, "ERROR: DeepTMHMM model weights under ${params.softwares_dir}/DeepTMHMM are incomplete. See message above."
    }
    def deeptmhmm_weights = file("${params.softwares_dir}/DeepTMHMM/DeepTMHMM-Academic-License-v1.0", checkIfExists: true)

    def ch_versions = channel.empty()

    //
    // Progress logging (nf-core/rnaseq-style): plain log.info lines announcing the
    // current stage and a running "done/total" count per sample. This exists
    // alongside -- not instead of -- Nextflow's own live process table, because that
    // table is a redrawing ANSI widget: it's illegible once written to a plain file,
    // which is exactly what happens on a `-bg`/detached run (see docs/usage.md). These
    // lines are just plain stdout, so they read fine in `.nextflow.log` / a redirected
    // log file, not only on an interactive terminal.
    //
    def pipeline_stages = [
        'FASTA_QC', 'DEEPCOIL2', 'PHOBIUS', 'INTERPROSCAN', 'INTERPROSCAN_MERGE',
        'DEEPLOC2', 'SIGNALP6', 'DEEPTMHMM', 'RGA_CLASSIFY', 'RGA_REPORT',
    ]
    def n_stages = pipeline_stages.size()
    // Same palette + params.monochrome_logs gate the nf-core completion summary already
    // uses (subworkflows/nf-core/utils_nfcore_pipeline's completionSummary) -- reused here
    // so progress logs match it and get disabled together with `--monochrome_logs`.
    def colors = logColours(params.monochrome_logs) as Map

    // Resolved asynchronously (channel.count() only emits once the upstream channel
    // closes) but in practice ready almost immediately: the samplesheet is a small,
    // already-in-memory list, not something fetched/streamed sample by sample.
    def total_samples = new java.util.concurrent.atomic.AtomicInteger(0)
    ch_samplesheet.count().subscribe { n ->
        total_samples.set(n as int)
        log.info ''
        log.info "${colors.bold}${colors.purple}${'=' * 60}${colors.reset}"
        log.info "${colors.bold}  RGAprofiler:${colors.reset} starting analysis of ${colors.bcyan}${n}${colors.reset} sample(s)"
        log.info "${colors.bold}  Stages${colors.reset} (${n_stages}): ${colors.dim}${pipeline_stages.join(' -> ')}${colors.reset}"
        log.info "${colors.bold}${colors.purple}${'=' * 60}${colors.reset}"
    }

    def stage_started  = new java.util.concurrent.ConcurrentHashMap()
    def stage_counters = new java.util.concurrent.ConcurrentHashMap()

    // Logs one "<sample> done (X/Y)" line per completed sample/chunk for `stage`, plus
    // a one-off "started" line the first time the stage is seen and a "complete" line
    // once every expected item has reported in. `total_override` lets INTERPROSCAN
    // (which fans out per FASTA chunk, not per sample) report against the chunk count
    // instead of the sample count; every other stage just uses total_samples.
    def logStageProgress = { String stage, String sample_id, Integer total_override = null ->
        def idx = pipeline_stages.indexOf(stage) + 1
        def tag = "${colors.bold}${colors.blue}[Stage ${idx}/${n_stages}]${colors.reset}"
        if (stage_started.putIfAbsent(stage, true) == null) {
            log.info ''
            log.info "${tag} ${colors.bold}${stage}${colors.reset} started"
        }
        def counter = stage_counters.computeIfAbsent(stage) { new java.util.concurrent.atomic.AtomicInteger(0) }
        def done = counter.incrementAndGet()
        def total = total_override != null ? total_override : total_samples.get()
        def total_display = total > 0 ? total.toString() : '?'
        def remaining = total > done ? total - done : 0
        log.info "${tag} ${stage}: ${colors.bcyan}${sample_id}${colors.reset} done (${colors.yellow}${done}/${total_display}${colors.reset} complete, ${colors.yellow}${remaining}${colors.reset} remaining)"
        if (total > 0 && done >= total) {
            log.info "${tag} ${colors.bold}${colors.green}${stage} complete (${done}/${total}) ✔${colors.reset}"
        }
    }

    //
    // Clean the input FASTA (dedup, strip stop-codon '*', uppercase) before handing
    // it to any prediction tool — real proteome FASTAs need this (e.g. trailing '*'
    // breaks Phobius/InterProScan). Chunk output feeds INTERPROSCAN's fan-out below;
    // no other tool consumes it (InterProScan alone gets a clear, measured benefit
    // from being split into many smaller jobs rather than one giant one).
    //
    FASTA_QC(ch_samplesheet, params.fasta_qc_chunk_size)
    ch_versions = ch_versions.mix(FASTA_QC.out.versions)
    FASTA_QC.out.fasta.subscribe { meta, fasta -> logStageProgress('FASTA_QC', meta.id) }

    //
    // Run per-tool prediction modules
    //
    DEEPCOIL2(FASTA_QC.out.fasta)
    ch_versions = ch_versions.mix(DEEPCOIL2.out.versions)
    DEEPCOIL2.out.results_dir.subscribe { meta, dir -> logStageProgress('DEEPCOIL2', meta.id) }

    PHOBIUS(FASTA_QC.out.fasta)
    ch_versions = ch_versions.mix(PHOBIUS.out.versions)
    PHOBIUS.out.predictions.subscribe { meta, preds -> logStageProgress('PHOBIUS', meta.id) }

    //
    // InterProScan scales much better as many smaller parallel jobs than one huge
    // single-FASTA job, so (unlike every other tool here) it runs once per FASTA_QC
    // chunk (.transpose() turns FASTA_QC's one-row-per-sample [meta, [chunk1,chunk2,...]]
    // into one row per chunk: [meta, chunk1], [meta, chunk2], ...), then INTERPROSCAN_MERGE
    // concatenates the per-chunk TSVs back into one file per sample (InterProScan's TSV
    // format has no header row, confirmed against a real run, so plain concatenation is
    // an exact, lossless merge) — nothing downstream needs to know chunking happened.
    //
    def ch_ipr_chunks = FASTA_QC.out.chunks.transpose()
    // Resolved asynchronously, same caveat as total_samples above: reads as '?' in the
    // progress log until every sample has finished FASTA_QC and the chunk count is final.
    def total_chunks = new java.util.concurrent.atomic.AtomicInteger(0)
    ch_ipr_chunks.count().subscribe { n -> total_chunks.set(n as int) }

    INTERPROSCAN(ch_ipr_chunks, ipr_dir)
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)
    INTERPROSCAN.out.tsv.subscribe { meta, tsv -> logStageProgress('INTERPROSCAN', meta.id, total_chunks.get()) }

    INTERPROSCAN_MERGE(INTERPROSCAN.out.tsv.groupTuple())
    INTERPROSCAN_MERGE.out.tsv.subscribe { meta, tsv -> logStageProgress('INTERPROSCAN_MERGE', meta.id) }

    DEEPLOC2(FASTA_QC.out.fasta, deeploc2_models, deeploc2_torch_cache)
    ch_versions = ch_versions.mix(DEEPLOC2.out.versions)
    DEEPLOC2.out.predictions.subscribe { meta, preds -> logStageProgress('DEEPLOC2', meta.id) }

    SIGNALP6(FASTA_QC.out.fasta, signalp6_models)
    ch_versions = ch_versions.mix(SIGNALP6.out.versions)
    SIGNALP6.out.predictions.subscribe { meta, preds -> logStageProgress('SIGNALP6', meta.id) }

    DEEPTMHMM(FASTA_QC.out.fasta, deeptmhmm_weights)
    ch_versions = ch_versions.mix(DEEPTMHMM.out.versions)
    DEEPTMHMM.out.predictions.subscribe { meta, preds -> logStageProgress('DEEPTMHMM', meta.id) }

    //
    // Classify proteins into RGA families/subclasses from the six tools' outputs above
    // (vendored from omatheuspimenta/rgapredictor, see docker/rga_classify/
    // and PLAN.md Stage 6). Joined by sample meta.id -- each upstream module emits exactly
    // one result per sample, so a plain join() is enough (no groupTuple needed).
    //
    def ch_deeptmhmm_gff3 = DEEPTMHMM.out.predictions
        .map { meta, files -> [ meta, files.find { it.toString().endsWith('.gff3') } ] }
    def ch_signalp6_predictions = SIGNALP6.out.predictions
        .map { meta, files -> [ meta, files.find { it.toString().contains('_predictions.txt') } ] }

    def ch_rga_classify_input = INTERPROSCAN_MERGE.out.tsv
        .join(PHOBIUS.out.predictions)
        .join(ch_deeptmhmm_gff3)
        .join(ch_signalp6_predictions)
        .join(DEEPLOC2.out.predictions)
        .join(DEEPCOIL2.out.results_dir)

    RGA_CLASSIFY(ch_rga_classify_input)
    ch_versions = ch_versions.mix(RGA_CLASSIFY.out.versions)
    RGA_CLASSIFY.out.results.subscribe { meta, files -> logStageProgress('RGA_CLASSIFY', meta.id) }

    //
    // Render the lightweight custom summary report (decision #3, PLAN.md -- no MultiQC).
    // RGA_REPORT reads a versions snapshot covering everything that ran *before* it
    // (FASTA_QC, the six tools, RGA_CLASSIFY); its own version is folded into the final
    // collated versions.yml below instead of being fed back into the report it generates,
    // which would be a circular dependency.
    //
    def ch_versions_for_report = softwareVersionsToYAML(ch_versions)
        .collectFile(sort: true, newLine: true)

    def ch_rga_report_input = RGA_CLASSIFY.out.results
        .map { meta, files ->
            [
                meta,
                files.find { it.toString().endsWith('rga_predictions.tsv') },
                files.find { it.toString().endsWith('rga_summary_counts.tsv') },
            ]
        }
        .combine(ch_versions_for_report)

    RGA_REPORT(ch_rga_report_input)
    ch_versions = ch_versions.mix(RGA_REPORT.out.versions)
    RGA_REPORT.out.report.subscribe { meta, files -> logStageProgress('RGA_REPORT', meta.id) }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'rgaprofiler_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
