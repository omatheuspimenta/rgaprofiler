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
include { DEEPCOIL2_MERGE    } from '../modules/local/deepcoil2_merge'
include { DEEPCOIL2_MERGE    as DEEPCOIL2_MERGE_BATCH    } from '../modules/local/deepcoil2_merge'
include { PHOBIUS            } from '../modules/local/phobius'
include { PHOBIUS_MERGE      } from '../modules/local/phobius_merge'
include { PHOBIUS_MERGE      as PHOBIUS_MERGE_BATCH      } from '../modules/local/phobius_merge'
include { INTERPROSCAN       } from '../modules/local/interproscan'
include { INTERPROSCAN_MERGE } from '../modules/local/interproscan_merge'
include { INTERPROSCAN_MERGE as INTERPROSCAN_MERGE_BATCH } from '../modules/local/interproscan_merge'
include { DEEPLOC2           } from '../modules/local/deeploc2'
include { DEEPLOC2_MERGE     } from '../modules/local/deeploc2_merge'
include { DEEPLOC2_MERGE     as DEEPLOC2_MERGE_BATCH     } from '../modules/local/deeploc2_merge'
include { SIGNALP6           } from '../modules/local/signalp6'
include { SIGNALP6_MERGE     } from '../modules/local/signalp6_merge'
include { SIGNALP6_MERGE     as SIGNALP6_MERGE_BATCH     } from '../modules/local/signalp6_merge'
include { DEEPTMHMM          } from '../modules/local/deeptmhmm'
include { DEEPTMHMM_MERGE    } from '../modules/local/deeptmhmm_merge'
include { DEEPTMHMM_MERGE    as DEEPTMHMM_MERGE_BATCH    } from '../modules/local/deeptmhmm_merge'
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
    if (params.use_gpu.toString() == 'auto') {
        def gpu_probe = ["bash", "${projectDir}/bin/detect_gpu.sh"].execute()
        signalp6_use_gpu = (gpu_probe.waitFor() == 0)
    } else {
        signalp6_use_gpu = (params.use_gpu.toString() == 'true')
    }
    // A forced GPU run on a host with no usable GPU fails every GPU task deep inside the
    // container runtime (docker: 'could not select device driver ""'). Only a warning, not
    // an error: under a grid executor the launch host legitimately has no GPU while the
    // compute nodes do.
    if (params.use_gpu.toString() == 'true' && ["bash", "${projectDir}/bin/detect_gpu.sh"].execute().waitFor() != 0) {
        log.warn "--use_gpu true (or -profile gpu) was set, but no usable NVIDIA GPU is visible on this host. " +
            "On the local executor every DeepCoil2/DeepLoc2/SignalP6/DeepTMHMM task will fail to start its container. " +
            "Use --use_gpu auto (default) or --use_gpu false to run them on CPU."
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
        'FASTA_QC', 'DEEPCOIL2', 'DEEPCOIL2_MERGE', 'PHOBIUS', 'PHOBIUS_MERGE', 'INTERPROSCAN', 'INTERPROSCAN_MERGE',
        'DEEPLOC2', 'DEEPLOC2_MERGE', 'SIGNALP6', 'SIGNALP6_MERGE', 'DEEPTMHMM', 'DEEPTMHMM_MERGE',
        'RGA_CLASSIFY', 'RGA_REPORT',
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
    // Clean the input FASTA (dedup, strip the trailing stop-codon '*', internal '*' -> X,
    // uppercase) before handing it to any prediction tool — real proteome FASTAs need
    // this (e.g. trailing '*' breaks Phobius/InterProScan) — and split it into sequence
    // blocks/chunks so DeepCoil2, Phobius, InterProScan, DeepLoc2, SignalP6 and DeepTMHMM
    // (below) can each fan out across them instead of ever processing an entire proteome
    // as one task.
    // --num_blocks (if set) takes priority and requests exactly that many chunks
    // (seqkit split2 --by-part balances sequences across them as evenly as possible);
    // otherwise chunks are sized by --fasta_qc_chunk_size sequences each, same as
    // before --num_blocks existed. Either way FASTA_QC always keeps whole FASTA
    // records intact -- chunking never splits a header from its sequence.
    // --num_blocks auto lets FASTA_QC pick that count itself from each sample's cleaned
    // sequence count and this host's CPUs (bin/fasta_qc.sh). Opt-in, not the default:
    // which proteins share a chunk slightly affects DeepCoil2 scores, SignalP6
    // probabilities and InterProScan's Gene3D boundaries (verified on real chunks), so
    // the chunking a run used is part of what makes its results reproducible.
    //
    def auto_blocks = params.num_blocks?.toString() == 'auto'
    def cpus_line = ["bash", "${projectDir}/bin/detect_host_resources.sh"].execute().text.readLines().find { line -> line.startsWith('cpus=') }
    def host_cpus = cpus_line ? (cpus_line.tokenize('=')[1] as int) : 1
    def split_mode  = auto_blocks ? 'auto' : (params.num_blocks ? 'parts' : 'size')
    def split_value = auto_blocks ? host_cpus : (params.num_blocks ?: params.fasta_qc_chunk_size)
    log.info ''
    log.info "${colors.bold}  Sequence batching:${colors.reset} " + (
        auto_blocks
            ? "--num_blocks ${colors.bcyan}auto${colors.reset} (chunk count chosen per sample from its sequence count and this host's ${host_cpus} CPUs; the exact count is in FASTA_QC's .command.err)"
            : params.num_blocks
                ? "--num_blocks ${colors.bcyan}${params.num_blocks}${colors.reset} (target chunk count, balanced by seqkit split2 --by-part)"
                : "--fasta_qc_chunk_size ${colors.bcyan}${params.fasta_qc_chunk_size}${colors.reset} sequences/chunk (--num_blocks not set)"
    )
    log.info "${colors.dim}  DeepCoil2, Phobius, InterProScan, DeepLoc2, SignalP6 and DeepTMHMM will each run once per resulting chunk.${colors.reset}"
    FASTA_QC(ch_samplesheet, split_mode, split_value)
    ch_versions = ch_versions.mix(FASTA_QC.out.versions)
    FASTA_QC.out.fasta.subscribe { meta, fasta -> logStageProgress.call('FASTA_QC', meta.id) }

    // .transpose() turns FASTA_QC's one-row-per-sample [meta, [chunk1,chunk2,...]]
    // into one row per chunk: [meta, chunk1], [meta, chunk2], ... -- shared by every
    // chunked tool below.
    def ch_chunks = FASTA_QC.out.chunks.transpose()
    // Resolved asynchronously, same caveat as total_samples above: reads as '?' in the
    // progress log until every sample has finished FASTA_QC and the chunk count is final.
    // Logged here (once) so a run's actual resulting chunk count -- not just the
    // requested --num_blocks/--fasta_qc_chunk_size value -- is visible in plain stdout,
    // confirming batching actually took effect (matters most for --num_blocks, whose
    // realised chunk count can differ slightly from the requested one -- see docs/output.md).
    def total_chunks = new java.util.concurrent.atomic.AtomicInteger(0)
    ch_chunks.count().subscribe { n ->
        total_chunks.set(n as int)
        log.info "${colors.bold}  Sequence batching:${colors.reset} produced ${colors.bcyan}${n}${colors.reset} chunk(s) total across ${colors.bcyan}${total_samples.get()}${colors.reset} sample(s)"
        // Advisory only -- chunking is never changed behind the user's back (it affects
        // some tools' numbers slightly, see above). InterProScan, the slowest per-chunk
        // step, uses 4 CPUs per chunk, so fewer chunks than 4-CPU slots leaves cores idle
        // on the machine running the pipeline (what matters for the local executor).
        def ips_slots = Math.max(1, host_cpus.intdiv(4))
        if (!params.num_blocks && n < ips_slots) {
            log.warn "Only ${n} chunk(s) for this ${host_cpus}-CPU machine: at most ${n} InterProScan task(s) (${n * 4} CPUs) can run at once. " +
                "--num_blocks auto would size chunks to use it fully (see docs/usage.md; the chunking a run uses slightly affects some tools' scores, so keep it fixed across runs you compare)."
        }
    }

    //
    // Run per-tool prediction modules. DeepCoil2, Phobius, InterProScan, DeepLoc2, SignalP6
    // and DeepTMHMM all run once per FASTA_QC chunk rather than once per sample -- these
    // six tools scale much better (and, for DeepCoil2 in particular, can fail outright
    // on a very large proteome) as many smaller parallel jobs than one huge single-FASTA
    // job. Each one's own *_MERGE module then reassembles the per-chunk outputs back
    // into one result per sample -- respecting each format's own record boundaries
    // (never an arbitrary line/byte split) -- so nothing downstream needs to know
    // chunking happened. (Phobius is single-threaded: unchunked, it took 8.8h on the
    // full R570 proteome, vs ~7min across 1000 chunks.)
    //
    // Every *_MERGE fan-in goes through splitForMerge (bottom of this file): a sample
    // with more chunks than one task can safely stage is first merged in batches
    // (*_MERGE_BATCH, same module), and those batch results are then merged again.
    //
    DEEPCOIL2(ch_chunks)
    ch_versions = ch_versions.mix(DEEPCOIL2.out.versions)
    DEEPCOIL2.out.results_dir.subscribe { meta, dir -> logStageProgress.call('DEEPCOIL2', meta.id, total_chunks.get()) }

    def deepcoil2_parts = splitForMerge(DEEPCOIL2.out.results_dir.groupTuple())
    DEEPCOIL2_MERGE_BATCH(deepcoil2_parts.batched)
    DEEPCOIL2_MERGE(deepcoil2_parts.direct.mix(DEEPCOIL2_MERGE_BATCH.out.results_dir.groupTuple()))
    DEEPCOIL2_MERGE.out.results_dir.subscribe { meta, dir -> logStageProgress.call('DEEPCOIL2_MERGE', meta.id) }

    PHOBIUS(ch_chunks)
    ch_versions = ch_versions.mix(PHOBIUS.out.versions)
    PHOBIUS.out.predictions.subscribe { meta, preds -> logStageProgress.call('PHOBIUS', meta.id, total_chunks.get()) }

    // Unlike the other merges, PHOBIUS_MERGE also takes the sample's clean FASTA: it
    // restores the record order an unchunked run writes (see modules/local/phobius_merge).
    // combine(by: 0), not join: a batched sample has one item per batch, all needing it.
    def phobius_parts = splitForMerge(PHOBIUS.out.predictions.groupTuple())
    PHOBIUS_MERGE_BATCH(phobius_parts.batched.combine(FASTA_QC.out.fasta, by: 0))
    PHOBIUS_MERGE(
        phobius_parts.direct
            .mix(PHOBIUS_MERGE_BATCH.out.predictions.groupTuple())
            .combine(FASTA_QC.out.fasta, by: 0)
    )
    PHOBIUS_MERGE.out.predictions.subscribe { meta, preds -> logStageProgress.call('PHOBIUS_MERGE', meta.id) }

    INTERPROSCAN(ch_chunks, ipr_dir)
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions)
    INTERPROSCAN.out.tsv.subscribe { meta, tsv -> logStageProgress.call('INTERPROSCAN', meta.id, total_chunks.get()) }

    def interproscan_parts = splitForMerge(INTERPROSCAN.out.tsv.groupTuple())
    INTERPROSCAN_MERGE_BATCH(interproscan_parts.batched)
    INTERPROSCAN_MERGE(interproscan_parts.direct.mix(INTERPROSCAN_MERGE_BATCH.out.tsv.groupTuple()))
    INTERPROSCAN_MERGE.out.tsv.subscribe { meta, tsv -> logStageProgress.call('INTERPROSCAN_MERGE', meta.id) }

    DEEPLOC2(ch_chunks, deeploc2_models, deeploc2_torch_cache)
    ch_versions = ch_versions.mix(DEEPLOC2.out.versions)
    DEEPLOC2.out.results_dir.subscribe { meta, dir -> logStageProgress.call('DEEPLOC2', meta.id, total_chunks.get()) }

    def deeploc2_parts = splitForMerge(DEEPLOC2.out.results_dir.groupTuple())
    DEEPLOC2_MERGE_BATCH(deeploc2_parts.batched)
    DEEPLOC2_MERGE(deeploc2_parts.direct.mix(DEEPLOC2_MERGE_BATCH.out.results_dir.groupTuple()))
    DEEPLOC2_MERGE.out.results_dir.subscribe { meta, dir -> logStageProgress.call('DEEPLOC2_MERGE', meta.id) }

    SIGNALP6(ch_chunks, signalp6_models)
    ch_versions = ch_versions.mix(SIGNALP6.out.versions)
    SIGNALP6.out.results_dir.subscribe { meta, dir -> logStageProgress.call('SIGNALP6', meta.id, total_chunks.get()) }

    def signalp6_parts = splitForMerge(SIGNALP6.out.results_dir.groupTuple())
    SIGNALP6_MERGE_BATCH(signalp6_parts.batched)
    SIGNALP6_MERGE(signalp6_parts.direct.mix(SIGNALP6_MERGE_BATCH.out.results_dir.groupTuple()))
    SIGNALP6_MERGE.out.results_dir.subscribe { meta, dir -> logStageProgress.call('SIGNALP6_MERGE', meta.id) }

    DEEPTMHMM(ch_chunks, deeptmhmm_weights)
    ch_versions = ch_versions.mix(DEEPTMHMM.out.versions)
    DEEPTMHMM.out.results_dir.subscribe { meta, dir -> logStageProgress.call('DEEPTMHMM', meta.id, total_chunks.get()) }

    def deeptmhmm_parts = splitForMerge(DEEPTMHMM.out.results_dir.groupTuple())
    DEEPTMHMM_MERGE_BATCH(deeptmhmm_parts.batched)
    DEEPTMHMM_MERGE(deeptmhmm_parts.direct.mix(DEEPTMHMM_MERGE_BATCH.out.results_dir.groupTuple()))
    DEEPTMHMM_MERGE.out.results_dir.subscribe { meta, dir -> logStageProgress.call('DEEPTMHMM_MERGE', meta.id) }

    //
    // Classify proteins into RGA families/subclasses from the six tools' outputs above
    // (vendored from omatheuspimenta/rgapredictor, see docker/rga_classify/
    // and PLAN.md Stage 6). Joined by sample meta.id -- each upstream module emits exactly
    // one result per sample, so a plain join() is enough (no groupTuple needed).
    //
    def ch_deeptmhmm_gff3 = DEEPTMHMM_MERGE.out.predictions
        .map { meta, files -> [ meta, files.find { it.toString().endsWith('.gff3') } ] }
    def ch_signalp6_predictions = SIGNALP6_MERGE.out.predictions
        .map { meta, files -> [ meta, files.find { it.toString().contains('_predictions.txt') } ] }

    // Chunking and CPU-vs-GPU both nudge some tools' numbers (docs/usage.md#chunking-and-reproducibility),
    // so RGA_CLASSIFY records them in run_metadata.json -- as requested *and* as resolved
    // (--num_blocks auto's real chunk count, --use_gpu auto's detected GPU), which
    // pipeline_info/params_*.json alone doesn't show. A single-chunk sample's chunks
    // output is a bare path rather than a list.
    def ch_run_settings = FASTA_QC.out.chunks
        .map { meta, chunks ->
            [
                meta,
                [
                    pipeline_version   : workflow.manifest.version,
                    num_blocks         : params.num_blocks,
                    fasta_qc_chunk_size: params.num_blocks ? null : params.fasta_qc_chunk_size,
                    n_chunks           : chunks instanceof List ? chunks.size() : 1,
                    use_gpu            : params.use_gpu,
                    gpu_used           : signalp6_use_gpu,
                ],
            ]
        }

    def ch_rga_classify_input = INTERPROSCAN_MERGE.out.tsv
        .join(PHOBIUS_MERGE.out.predictions)
        .join(ch_deeptmhmm_gff3)
        .join(ch_signalp6_predictions)
        .join(DEEPLOC2_MERGE.out.predictions)
        .join(DEEPCOIL2_MERGE.out.results_dir)
        .join(ch_run_settings)

    RGA_CLASSIFY(ch_rga_classify_input)
    ch_versions = ch_versions.mix(RGA_CLASSIFY.out.versions)
    RGA_CLASSIFY.out.results.subscribe { meta, files -> logStageProgress.call('RGA_CLASSIFY', meta.id) }

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
    RGA_REPORT.out.report.subscribe { meta, files -> logStageProgress.call('RGA_REPORT', meta.id) }

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Routes one tool's per-sample [meta, [chunk outputs]] fan-in to its *_MERGE step.
// A single task can't stage an unbounded number of inputs: Nextflow writes one
// 'rm -f'/'ln -s' pair per input into a single bash function (nxf_stage) in
// .command.run, and bash itself segfaults (exit 139) parsing that function once it
// holds somewhere past ~10k inputs -- before the task's own script ever runs, and
// with no executor setting to avoid it. So:
//   - direct:  samples with <= batch_size chunks go straight to the final *_MERGE,
//              exactly as before.
//   - batched: larger samples are cut into [meta, batch] items of <= batch_size
//              chunks for *_MERGE_BATCH; its per-batch results are then grouped
//              again and merged by the final *_MERGE (every *_MERGE accepts its own
//              output shape as input). batch_size^2 = 1,000,000 chunks per sample
//              before even the final merge would exceed batch_size inputs.
// Items are sorted by path so batch membership -- and so every merge task's input
// hash -- is identical across runs, keeping -resume effective for these steps.
//
def splitForMerge(ch_grouped) {
    def batch_size = 1000
    return ch_grouped
        .flatMap { meta, items ->
            def sorted = items.sort(false) { item -> item.toString() }
            sorted.size() <= batch_size
                ? [[meta, sorted, 'direct']]
                : sorted.collate(batch_size).collect { batch -> [meta, batch, 'batched'] }
        }
        .branch { meta, items, route ->
            direct: route == 'direct'
                return [meta, items]
            batched: true
                return [meta, items]
        }
}
