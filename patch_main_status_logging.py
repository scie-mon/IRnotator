#!/usr/bin/env python3
from pathlib import Path

src = Path('main.nf')
dst = Path('main.status.nf')
text = src.read_text()

old_init = '''    initial_deeptmhmm_res = INITIALIZE_DEEPTMHMM_RESULTS(
        hmm_hit_faa_ch,
        Channel.value(salvage_paths)
    )

    def generated_results_ch
'''
new_init = '''    initial_deeptmhmm_res = INITIALIZE_DEEPTMHMM_RESULTS(
        hmm_hit_faa_ch,
        Channel.value(salvage_paths)
    )

    initial_manifest_ch = initial_deeptmhmm_res.manifest.map { manifest_tsv ->
        def rows = manifest_tsv.toFile().readLines().drop(1).findAll { it }.collect {
            it.split('\\t', -1)
        }
        def salvaged = rows.count { fields -> fields[3] == 'salvaged' }
        log.info "DeepTMHMM candidates: ${rows.size()}; salvaged: ${salvaged}; queued for new prediction: ${rows.size() - salvaged}"
        manifest_tsv
    }

    def generated_results_ch
'''
old_final = '''    final_deeptmhmm_res = FINALIZE_DEEPTMHMM_RESULTS(
        initial_deeptmhmm_res.manifest,
        generated_results_ch
    )
'''
new_final = '''    final_deeptmhmm_res = FINALIZE_DEEPTMHMM_RESULTS(
        initial_manifest_ch,
        generated_results_ch
    )
'''

if text.count(old_init) != 1 or text.count(old_final) != 1:
    raise SystemExit('Expected main.nf blocks were not found exactly once; no file was written.')

dst.write_text(text.replace(old_init, new_init).replace(old_final, new_final))
