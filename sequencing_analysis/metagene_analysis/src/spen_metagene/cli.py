from __future__ import annotations

import argparse, json
from pathlib import Path
import yaml

from .annotation import save_models, save_gene_covariates
from .classes import classify_file, early_agreement
from .clap import count_sample, combine_sample_tables
from .counts import prepare_counts, merge_de_annotation, write_manifest
from .historical import summarize
from .historical_profiles import build as build_historical_profiles, compare_stability
from .nascent_counter import count as count_nascent
from .plotting import paired_volcano, ma_plot
from .report import build_report
from .suite import make_suite
from .summaries import build as build_summaries
from .validate import validate_config
from .provenance import write_versions
from .qc import build as build_qc
from .audit import write_audit
from .distributions import plot_distribution_prototype
from .matched import run as run_matched
from .matched_report import build as build_matched_report
from .expression_tertiles import run as run_expression_tertiles
from .artifact_validation import run as run_artifact_validation


def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest="cmd",required=True)
    q=sub.add_parser("count-nascent"); q.add_argument("--bam"); q.add_argument("--gtf"); q.add_argument("--mask"); q.add_argument("--label"); q.add_argument("--output")
    q=sub.add_parser("prepare-counts"); q.add_argument("--expr-dir"); q.add_argument("--late-dir"); q.add_argument("--output-dir")
    q=sub.add_parser("merge-annotation"); q.add_argument("--de"); q.add_argument("--annotation"); q.add_argument("--output")
    q=sub.add_parser("classify"); q.add_argument("--input"); q.add_argument("--output"); q.add_argument("--delta",type=float,default=0.2630344058337938)
    q=sub.add_parser("agreement"); q.add_argument("--one"); q.add_argument("--two"); q.add_argument("--output-prefix"); q.add_argument("--delta",type=float,default=.2630344058337938)
    q=sub.add_parser("volcano"); q.add_argument("--input"); q.add_argument("--output"); q.add_argument("--class-col"); q.add_argument("--title"); q.add_argument("--delta",type=float); q.add_argument("--fdr",type=float)
    q=sub.add_parser("ma"); q.add_argument("--input"); q.add_argument("--output"); q.add_argument("--class-col"); q.add_argument("--title")
    q=sub.add_parser("models"); q.add_argument("--gtf"); q.add_argument("--output")
    q=sub.add_parser("gene-covariates"); q.add_argument("--gtf"); q.add_argument("--output")
    q=sub.add_parser("count-clap"); q.add_argument("--models"); q.add_argument("--bam"); q.add_argument("--sample"); q.add_argument("--output-prefix"); q.add_argument("--min-mapq",type=int,default=20)
    q=sub.add_parser("combine-clap"); q.add_argument("--inputs",nargs="+"); q.add_argument("--samples",nargs="+"); q.add_argument("--output")
    q=sub.add_parser("manifest"); q.add_argument("--config"); q.add_argument("--output")
    q=sub.add_parser("validate"); q.add_argument("--config"); q.add_argument("--output")
    q=sub.add_parser("historical"); q.add_argument("--input-dir"); q.add_argument("--output")
    q=sub.add_parser("historical-profiles"); q.add_argument("--input-dir"); q.add_argument("--models"); q.add_argument("--output-dir")
    q=sub.add_parser("ebs-stability"); q.add_argument("--archived-dir"); q.add_argument("--stability-dir"); q.add_argument("--permutations",type=int); q.add_argument("--output")
    q=sub.add_parser("suite"); q.add_argument("--signals-json"); q.add_argument("--classes-json"); q.add_argument("--agreements-json"); q.add_argument("--covariates"); q.add_argument("--models"); q.add_argument("--output-dir"); q.add_argument("--bootstrap",type=int,default=100); q.add_argument("--seed",type=int,default=20260812)
    q=sub.add_parser("report"); q.add_argument("--output-root"); q.add_argument("--output"); q.add_argument("--trimmed-output")
    q=sub.add_parser("summaries"); q.add_argument("--classes-json"); q.add_argument("--output-dir")
    q=sub.add_parser("versions"); q.add_argument("--output")
    q=sub.add_parser("qc"); q.add_argument("--root")
    q=sub.add_parser("audit"); q.add_argument("--root"); q.add_argument("--output")
    q=sub.add_parser("distribution-prototype"); q.add_argument("--features"); q.add_argument("--classes"); q.add_argument("--output-prefix"); q.add_argument("--title",default="Aux 2 h vs control — total RNA"); q.add_argument("--min-kde-n",type=int,default=30)
    q=sub.add_parser("expression-match"); q.add_argument("--features"); q.add_argument("--classes"); q.add_argument("--output-prefix"); q.add_argument("--title"); q.add_argument("--bin-width",type=float,default=.1); q.add_argument("--seed",type=int,default=20260812); q.add_argument("--min-kde-n",type=int,default=30)
    q=sub.add_parser("matched-report"); q.add_argument("--output-root"); q.add_argument("--output")
    q=sub.add_parser("expression-tertiles"); q.add_argument("--output-root"); q.add_argument("--output"); q.add_argument("--seed",type=int,default=20260812); q.add_argument("--bootstrap",type=int,default=100); q.add_argument("--min-kde-n",type=int,default=30)
    q=sub.add_parser("artifact-validation"); q.add_argument("--output-root"); q.add_argument("--config",required=True)
    a=p.parse_args()
    if a.cmd=="count-nascent": count_nascent(a.bam,a.gtf,a.mask,a.label,a.output)
    elif a.cmd=="prepare-counts": prepare_counts(a.expr_dir,a.late_dir,a.output_dir)
    elif a.cmd=="merge-annotation": merge_de_annotation(a.de,a.annotation,a.output)
    elif a.cmd=="classify": classify_file(a.input,a.output,a.delta)
    elif a.cmd=="agreement":
        import pandas as pd
        one=pd.read_csv(a.one,sep="\t"); two=pd.read_csv(a.two,sep="\t")
        for f in [.05,.10]: early_agreement(one,two,f,a.delta).to_csv(f"{a.output_prefix}.fdr{int(f*100):02d}.tsv",sep="\t",index=False)
    elif a.cmd=="volcano": paired_volcano(a.input,a.output,a.class_col,a.title,a.delta,a.fdr)
    elif a.cmd=="ma": ma_plot(a.input,a.output,a.class_col,a.title)
    elif a.cmd=="models": save_models(a.gtf,a.output)
    elif a.cmd=="gene-covariates": save_gene_covariates(a.gtf,a.output)
    elif a.cmd=="count-clap": count_sample(a.models,a.bam,a.sample,a.output_prefix,a.min_mapq)
    elif a.cmd=="combine-clap": combine_sample_tables(a.inputs,a.samples,a.output)
    elif a.cmd=="manifest":
        with open(a.config) as fh: cfg=yaml.safe_load(fh)
        write_manifest(cfg,a.output)
    elif a.cmd=="validate":
        with open(a.config) as fh: cfg=yaml.safe_load(fh)
        validate_config(cfg,a.output)
    elif a.cmd=="historical": summarize(a.input_dir,a.output)
    elif a.cmd=="historical-profiles": build_historical_profiles(a.input_dir,a.models,a.output_dir)
    elif a.cmd=="ebs-stability": compare_stability(a.archived_dir,a.stability_dir,a.permutations,a.output)
    elif a.cmd=="suite":
        with open(a.signals_json) as fh: signals=json.load(fh)
        with open(a.classes_json) as fh: classes=json.load(fh)
        agreements={}
        if a.agreements_json:
            with open(a.agreements_json) as fh: agreements=json.load(fh)
        make_suite(signals,classes,a.output_dir,a.bootstrap,a.seed,agreements,a.covariates,a.models)
    elif a.cmd=="report": build_report(a.output_root,a.output,a.trimmed_output)
    elif a.cmd=="summaries":
        with open(a.classes_json) as fh: files=json.load(fh)
        build_summaries(files,a.output_dir)
    elif a.cmd=="versions": write_versions(a.output)
    elif a.cmd=="qc": build_qc(a.root)
    elif a.cmd=="audit": write_audit(a.root,a.output)
    elif a.cmd=="distribution-prototype": plot_distribution_prototype(a.features,a.classes,a.output_prefix,a.title,a.min_kde_n)
    elif a.cmd=="expression-match": run_matched(a.features,a.classes,a.output_prefix,a.title,a.bin_width,a.seed,a.min_kde_n)
    elif a.cmd=="matched-report": build_matched_report(a.output_root,a.output)
    elif a.cmd=="expression-tertiles": run_expression_tertiles(a.output_root,a.output,a.seed,a.bootstrap,a.min_kde_n)
    elif a.cmd=="artifact-validation": run_artifact_validation(a.output_root,a.config)


if __name__=="__main__": main()
