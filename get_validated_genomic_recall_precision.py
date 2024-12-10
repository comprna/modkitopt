import argparse
import pickle
import time

from types import SimpleNamespace


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_bed",
                        help="Input bed file, lifted to genome with columns:\n"
                            "chr pos0 pos1    .... coverage stoichiometry site_confidence",
                        required=True)
    parser.add_argument("-v", "--validated",
                        help="Input pickle file with a set of validated sites:\n"
                            "{chr_pos}",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="Input bed file with validated sites:\n"
                            "chr pos0 pos1 ",
                        required=True)
    parser.add_argument("-d", "--discard",
                        help='Path to set with sites to discard from analysis',
                        type=int,
                        default=None,
                        required=False)
    parser.add_argument("-a", "--aggregate",
                        choices=["max", "avg"],
                        help='Method for handling aggregated transcripts at same genomic location.\n'
                            '<max/avg> ; Default is max',
                        default= "max",
                        required=False)
    parser.add_argument("-m", "--metric",
                        choices=["model2", "rate"],
                        help='Metric for selecting sites.'
                            '<model2/rate> ; Default is model2',
                        default= "model2",
                        required=False)
    parser.add_argument("--min_stoich",
                        help='Minimum stoichiometry for selecting sites. (Between 0 and 1)',
                        default= -1.0,
                        type=float,
                        required=False)
    # args = parser.parse_args()

    # Local testing
    top_dir = "/home/alex/OneDrive/Projects/m6A_proteins/1_prelim_analysis/1_preprocess_m6A/drs/1_optimise_modkit"
    args = SimpleNamespace()
    args.validated = "./example/validated_genomic_sites.pickle"
    args.input_bed = "./example/predicted_genomic_sites.bed"
    args.output = "./example/test_output.tsv"
    args.aggregate = "avg"
    args.metric = "rate"
    args.discard = None
    args.min_stoich = 0.8

    st = time.time()

    # Load list of sites to discard
    if args.discard:
        with open(args.discard,"rb") as p:
            discard_sites = pickle.load(p)
    else:
        discard_sites = set()

    # Load list of validated sites
    with open(args.validated,"rb") as p:
        validated_sites = pickle.load(p)

    # Parse site level info to and store site thesholds
    preds_dct_val = {}
    preds_dct_all = {}
    with open(args.input_bed) as f:
        f.readline() # Pass over header
        for line in f:
            fields = line.strip().split("\t")
            coverage, stoich, prob = fields[-3:]
            contig, start, end = fields[:3]
            site_index = f"{contig}_{end}"

            # If sites are selected based on modification rate, then... TODO:stefan what is logic here?
            if args.metric == "rate":
                prob = stoich

            # Some m6A stoichiometry predictors output None if... TODO:stefan
            try:
                if float(stoich) > args.min_stoich:
                    p_predicted = float(prob)
                else:
                    p_predicted = 0
            except:
                print(stoich, prob, "skipped")
                continue

            # Skip over this site if it needs to be discarded
            if site_index in discard_sites:
                continue

            # Aggregate transcriptomic sites corresponding to the same genomic
            # site by taking the transcriptomic site with maximum stoichiometry
            if args.aggregate == "max":
                if site_index in preds_dct_all:
                    preds_dct_all[site_index] = max(preds_dct_all[site_index], p_predicted)
                else:
                    preds_dct_all[site_index] = p_predicted
                if site_index in validated_sites:
                    if site_index in preds_dct_val:
                            preds_dct_val[site_index] = max(preds_dct_val[site_index], p_predicted)
                    else:
                        preds_dct_val[site_index] = p_predicted

            # Aggregate transcriptomic sites corresponding to the same genomic
            # site by taking the average stoichiometry across transcriptomic
            # sites, accounting for transcript coverage
            if args.aggregate == "avg":
                coverage = float(coverage)
                if site_index in preds_dct_all:
                    predsum, covsum = preds_dct_all[site_index]
                    predsum += p_predicted*coverage
                    covsum += coverage
                    preds_dct_all[site_index] = [predsum, covsum] # store avg weighted probability if multiple transcripts
                else:
                    preds_dct_all[site_index] = [p_predicted*coverage, coverage]

                if site_index in validated_sites:
                    if site_index in preds_dct_val:
                        predsum, covsum = preds_dct_val[site_index]
                        predsum += p_predicted * coverage
                        covsum += coverage
                        preds_dct_val[site_index] = [predsum, covsum]  # store avg weighted probability if multiple transcripts
                    else:
                        preds_dct_val[site_index] = [p_predicted*coverage, coverage]

    # Get list of site probabilities
    if args.aggregate == "max":
        p_lst_all = [item for key, item in preds_dct_all.items()]
        p_lst_val = [item for key, item in preds_dct_val.items()]
    if args.aggregate == "avg":
        p_lst_all = [item[0] / item[1] for key, item in preds_dct_all.items()]
        p_lst_val = [item[0] / item[1] for key, item in preds_dct_val.items()]

    # Sort list of site probabilities
    p_lst_val = sorted(p_lst_val)
    p_lst_all = sorted(p_lst_all)

    # Initialise thresholds
    thresholds = [0, 1/10**7, 1/10**6, 1/10**5, 1/10**4] + [x/1000 for x in range(1, 900)]
    
    # Add range of thresholds
    # 0.9, 0.901, 0.902 ... 0.989, 0.990, 0.9901 ... 0.9990 ... 0.99999999989 ... 1
    base_n_lst = [0.9]
    for n1 in range(1, 10):
        base_n  = str(base_n_lst[-1])
        for n2 in range(0, 9):
            base_n_to_add = float(base_n + str(n2)) # TODO: String concatenation in this manner is slow
            base_n_lst.append(base_n_to_add)
            for n3 in range(0, 10):
                threshold_to_add = float(base_n + str(n2) + str(n3))
                thresholds.append(threshold_to_add)
        base_n_lst.append(float(base_n + str(9)))
    thresholds.append(1)

    # Get n predicted at thresholds for validated sites
    # TODO: Just compares number predicted, not whether each individual site is correctly predicted
    index_p, index_t = 0, 0
    out_vals = []
    while index_p < len(p_lst_val) and index_t < len(thresholds):
        prob = p_lst_val[index_p]
        thresh = thresholds[index_t]
        if prob < thresh:
            index_p += 1
        else:
            index_t += 1
            out_vals.append([thresh, 1 - index_p/len(p_lst_val), len(p_lst_val) - index_p])

    while index_t < len(thresholds):
        thresh = thresholds[index_t]
        out_vals.append([thresh, 0, 0])
        index_t += 1

    # Get n predicted at thresholds for all sites
    index_p, index_t = 0, 0
    out_vals_all = []
    while index_p < len(p_lst_all) and index_t < len(thresholds):
        prob = p_lst_all[index_p]
        thresh = thresholds[index_t]
        if prob < thresh:
            index_p += 1
        else:
            index_t += 1
            out_vals_all.append([thresh, 1 - index_p/len(p_lst_all), len(p_lst_all) - index_p])

    while index_t < len(thresholds):
        thresh = thresholds[index_t]
        out_vals_all.append([thresh, 0, 0])
        index_t += 1

    # Write metrics to file
    with open(args.output, "w+") as of:
        of.write("site_threshold\tvalidated_called\tvalidated_rate\tall_called\tall_rate\tvalidated_precision\n")
        for i, row in enumerate(out_vals):
            thresh, validated_freq, validated_pred = row
            thresh, all_freq, all_pred = out_vals_all[i]
            if all_pred > 0:
                called_validated_percent = validated_pred / all_pred
            else:
                called_validated_percent = 0
            of.write("\t".join([str(x) for x in [thresh, validated_pred, validated_freq, all_pred, all_freq, called_validated_percent]]) + "\n")

    # Note the time taken to run this analysis
    print("all done in ", time.time() - st, "seconds") # TODO: Convert to f-string


if __name__ == "__main__":
    main()
