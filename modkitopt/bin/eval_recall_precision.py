import argparse
import pickle


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
            seqname = fields[0]
            site = fields[1]
            truth_sites.add(f"{seqname}_{site}")

    # Parse predicted sites in input BED file
    preds_validated = {} # Predicted sites that are also in the ground truth sites
    preds_all = {} # All predicted sites
    with open(args.input_bed) as f:
        f.readline() # Pass over header
        for line in f:
            fields = line.strip().split("\t")
            coverage, stoich = fields[-3:-1] # TODO: [-2:]
            chromosome, start, end = fields[:3]
            pred_site = f"chr{chromosome}_{end}"

            # Some m6A stoichiometry predictors output None
            try:
                stoich = float(stoich)
            except:
                print(f"Skipped: {line}")
                continue

            # Aggregate transcriptomic sites corresponding to the same genomic
            # site by taking the average stoichiometry across transcriptomic
            # sites, weighted by transcript coverage
            coverage = float(coverage)
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
    thresholds = [0, 1/10**7, 1/10**6, 1/10**5, 1/10**4] + [x/1000 for x in range(1, 900)]
    
    # Add range of thresholds
    # 0.9, 0.901, 0.902 ... 0.989, 0.990, 0.9901 ... 0.9990 ... 0.99999999989 ... 1
    base_n_lst = [0.9]
    for n1 in range(1, 10):
        base_n  = str(base_n_lst[-1])
        for n2 in range(0, 9):
            base_n_to_add = float(base_n + str(n2))
            base_n_lst.append(base_n_to_add)
            for n3 in range(0, 10):
                threshold_to_add = float(base_n + str(n2) + str(n3))
                thresholds.append(threshold_to_add)
        base_n_lst.append(float(base_n + str(9)))
    thresholds.append(1)

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
            out_vals.append([thresh, 1 - index_p/len(p_validated), len(p_validated) - index_p])

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
            out_vals_all.append([thresh, 1 - index_p/len(p_all), len(p_all) - index_p])

    while index_t < len(thresholds):
        thresh = thresholds[index_t]
        out_vals_all.append([thresh, 0, 0])
        index_t += 1

    # Write metrics to file
    with open(args.output, "w+") as of:
        # threshold   --> stoichiometry threshold to consider a site as modified
        # recall      --> proportion of validated sites that were predicted
        # precision   --> proportion of predicted sites that were validated
        of.write("threshold\trecall\tprecision\n")
        for i, row in enumerate(out_vals):
            thresh, recall, tp = row
            thresh, _, all_pred = out_vals_all[i]
            if all_pred > 0:
                precision = tp / all_pred
            else:
                precision = 0
            of.write("\t".join([str(x) for x in [thresh, recall, precision]]) + "\n")


if __name__ == "__main__":
    main()
