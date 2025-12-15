#!/usr/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_bed",
                        help="BED file of predicted sites (genomic) with columns:\n"
                            "chr pos0 pos1 .... coverage stoichiometry site_confidence",
                        required=True)
    parser.add_argument("-t", "--truth",
                        help="Ground truth sites (genomic) with first two columns:\n"
                            "chr pos",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="Predicted sites validated by ground truth, with columns:\n"
                            "chr pos0 pos1",
                        required=True)
    args = parser.parse_args()

    # Load ground truth sites
    truth_sites = set()
    with open(args.truth,"r") as f:
        f.readline()
        for line in f:
            fields = line.split("\t")
            seqname = fields[0].lstrip("chr")
            site = fields[1]
            truth_sites.add(f"{seqname}_{site}")

    # Parse predicted sites in input BED file
    preds_validated = {} # Predicted sites that are also in the ground truth sites
    preds_all = {} # All predicted sites
    with open(args.input_bed) as f:
        f.readline() # Pass over header
        for line in f:
            fields = line.strip().split("\t")
            coverage, stoich = fields[-2:]
            chromosome, start, end = fields[:3]
            chromosome = chromosome.lstrip("chr")
            pred_site = f"{chromosome}_{end}"

            # Only consider sites with coverage at least 20 reads
            coverage = int(coverage)
            if coverage < 20:
                continue

            # Some m6A stoichiometry predictors output None
            try:
                stoich = float(stoich)
            except:
                print(f"Skipped: {line}")
                continue

            # Convert stoichiometry from % to fraction
            stoich = stoich / 100

            # Aggregate transcriptomic sites corresponding to the same genomic
            # site by taking the average stoichiometry across transcriptomic
            # sites, weighted by transcript coverage
            if pred_site in preds_all:
                pred_sum, cov_sum = preds_all[pred_site]
                pred_sum += stoich * coverage
                cov_sum += coverage
                preds_all[pred_site] = [pred_sum, cov_sum]
            else:
                preds_all[pred_site] = [stoich * coverage, coverage]

            if pred_site in truth_sites:
                if pred_site in preds_validated:
                    pred_sum, cov_sum = preds_validated[pred_site]
                    pred_sum += stoich * coverage
                    cov_sum += coverage
                    preds_validated[pred_site] = [pred_sum, cov_sum]
                else:
                    preds_validated[pred_site] = [stoich * coverage, coverage]

    # Get list of site probabilities
    p_all = [item[0] / item[1] for key, item in preds_all.items()]
    p_validated = [item[0] / item[1] for key, item in preds_validated.items()]

    # Sort list of site probabilities
    p_validated = sorted(p_validated)
    p_all = sorted(p_all)

    # Initialise stoichiometry thresholds
    thresholds = [t / 1000 for t in range(0, 1000)]

    # Get n validated sites predicted per threshold
    index_p, index_t = 0, 0
    out_vals = []
    while index_p < len(p_validated) and index_t < len(thresholds):
        prob = p_validated[index_p]
        thresh = thresholds[index_t]
        if prob < thresh:
            index_p += 1
        else:
            index_t += 1
            n_true_positives = len(p_validated) - index_p # No. of predictions that exceed the threshold
            recall = n_true_positives / len(truth_sites)
            out_vals.append([thresh, recall, n_true_positives])

    while index_t < len(thresholds):
        thresh = thresholds[index_t]
        out_vals.append([thresh, 0, 0])
        index_t += 1

    # Get n sites (total, not just validated sites) predicted per threshold
    index_p, index_t = 0, 0
    out_vals_all = []
    while index_p < len(p_all) and index_t < len(thresholds):
        prob = p_all[index_p]
        thresh = thresholds[index_t]
        if prob < thresh:
            index_p += 1
        else:
            index_t += 1
            n_predicted_positives = len(p_all) - index_p
            out_vals_all.append([thresh, n_predicted_positives])

    while index_t < len(thresholds):
        thresh = thresholds[index_t]
        out_vals_all.append([thresh, 0])
        index_t += 1

    # Write metrics to file
    with open(args.output, "w+") as out:
        # threshold   --> stoichiometry threshold to consider a site as modified
        # recall      --> proportion of validated sites that were predicted
        # precision   --> proportion of predicted sites that were validated
        out.write("threshold\trecall\tprecision\n")
        for i, row in enumerate(out_vals):
            thresh, recall, n_true_positives = row
            thresh, n_predicted_positives = out_vals_all[i]
            if n_predicted_positives > 0:
                precision = n_true_positives / n_predicted_positives
            else:
                precision = 0
            out.write(f"{thresh}\t{recall:.3f}\t{precision:.3f}\n")


if __name__ == "__main__":
    main()
