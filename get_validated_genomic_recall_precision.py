import cProfile
import time
import pickle
import argparse

from types import SimpleNamespace


def write_output(out_path, out_vals_val, out_vals_all):
    with open(out_path, "w+") as of:
        of.write("site_threshold\tvalidated_called\tvalidated_rate\tall_called\tall_rate\tvalidated_precision\n")
        for i, row in enumerate(out_vals_val):
            thresh, validated_freq, validated_pred = row
            thresh, all_freq, all_pred = out_vals_all[i]
            if all_pred > 0:
                called_validated_percent = validated_pred / all_pred
            else:
                called_validated_percent = 0
            of.write("\t".join([str(x) for x in [thresh, validated_pred, validated_freq, all_pred, all_freq, called_validated_percent]]) + "\n")


def get_n_predicted(p_lst, thresholds):
    # two pointers to get n predicted at thresholds for validated sites
    # TODO: Just compares number predicted, not whether each individual site is correctly predicted
    index_p, index_t = 0, 0
    out_vals = []
    while index_p < len(p_lst) and index_t < len(thresholds):
        prob = p_lst[index_p]
        thresh = thresholds[index_t]
        if prob < thresh:
            index_p += 1
        else:
            index_t += 1
            out_vals.append([thresh, 1 - index_p/len(p_lst), len(p_lst) - index_p])

    while index_t < len(thresholds):
        thresh = thresholds[index_t]
        out_vals.append([thresh, 0, 0])
        index_t += 1

    return out_vals


def define_thresholds():
    # Initialise thresholds
    thresholds = [0, 1/10**7, 1/10**6, 1/10**5, 1/10**4] + [x/1000 for x in range(1, 900)]
    # Add range of thresholds 0.9, 0.901, 0.902 ... 0.989, 0.990, 0.9901 ... 0.9990 ... 0.99999999989 ... 1
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
    return thresholds


def get_probs_list(preds_dct, aggregate):
    if aggregate == "max":
        p_lst = [item for key, item in preds_dct.items()]

    if aggregate == "avg":
        p_lst = [item[0] / item[1] for key, item in preds_dct.items()]

    return sorted(p_lst)


def parse_sites(predicted_path, validated_set, metric, min_stoich, aggregate, discard_set):
    preds_dct_val = {}
    preds_dct_all = {}
    with open(predicted_path) as f:
        f.readline()
        for line in f:
            line_lst = line.strip().split("\t")
            coverage, stoich, prob = line_lst[-3:]
            contig, start, end = line_lst[:3]
            site_index = f"{contig}_{end}"

            if metric == "rate":
                prob = stoich

            try:
                if float(stoich) > min_stoich:
                    p_predicted = float(prob)
                else:
                    p_predicted = 0
            except:
                print(stoich, prob, "skipped")
                continue

            if site_index in discard_set:
                continue

            if aggregate == "max":
                if site_index in preds_dct_all:
                    preds_dct_all[site_index] = max(preds_dct_all[site_index], p_predicted) # store maximum probability if multiple transcripts
                else:
                    preds_dct_all[site_index] = p_predicted # store maximum probability if multiple transcripts
                if site_index in validated_set:
                    if site_index in preds_dct_val:    # store maximum probability if multiple transcripts
                            preds_dct_val[site_index] = max(preds_dct_val[site_index], p_predicted)
                    else:
                        preds_dct_val[site_index] = p_predicted

            if aggregate == "avg":
                coverage = float(coverage)
                if site_index in preds_dct_all:
                    predsum, covsum = preds_dct_all[site_index]
                    predsum += p_predicted*coverage
                    covsum += coverage
                    preds_dct_all[site_index] = [predsum, covsum] # store avg weighted probability if multiple transcripts
                else:
                    preds_dct_all[site_index] = [p_predicted*coverage, coverage]

                if site_index in validated_set:
                    if site_index in preds_dct_val:
                        predsum, covsum = preds_dct_val[site_index]
                        predsum += p_predicted * coverage
                        covsum += coverage
                        preds_dct_val[site_index] = [predsum, covsum]  # store avg weighted probability if multiple transcripts
                    else:
                        preds_dct_val[site_index] = [p_predicted*coverage, coverage]
    return preds_dct_all, preds_dct_val


def main():
    parser = argparse.ArgumentParser()
    OPTIONAL = parser._action_groups.pop()
    REQUIRED = parser.add_argument_group('required arguments')

    REQUIRED.add_argument("-i", "--input_bed",
                        help="Input bed file, lifted to genome with columns:\n"
                            "chr pos0 pos1    .... coverage stoichiometry site_confidence",
                        required=True)

    REQUIRED.add_argument("-v", "--validated",
                        help="Input pickle file with a set of validated sites:\n"
                            "{chr_pos}",
                        required=True)

    REQUIRED.add_argument("-o", "--output",
                        help="Input bed file with validated sites:\n"
                            "chr pos0 pos1 ",
                        required=True)

    OPTIONAL.add_argument("-d", "--discard",
                        help='Path to set with sites to discard from analysis',
                        type=int,
                        default=None
                        )

    OPTIONAL.add_argument("-a", "--aggregate",
                        help='Method for handling aggregated transcripts at same genomic location.\n'
                            '<max/avg> ; Default is max',
                        default= "max"
                        )

    OPTIONAL.add_argument("-m", "--metric",
                        help='Metric for selecting sites.'
                            '<model2/rate> ; Default is model2',
                        default= "model2"
                        )

    OPTIONAL.add_argument("--min_stoich",
                        help='Minimum stoichiometry for selecting sites. (Between 0 and 1)',
                        default= -1.0,
                        type=float
                        )

    parser._action_groups.append(OPTIONAL)

    # ARGS = parser.parse_args()
    top_dir = "/home/alex/OneDrive/Projects/m6A_proteins/1_prelim_analysis/1_preprocess_m6A/drs/1_optimise_modkit"
    ARGS = SimpleNamespace()
    # ARGS.validated = "./example/validated_genomic_sites.pickle"
    # ARGS.input_bed = "./example/predicted_genomic_sites.bed"
    # ARGS.output = "./example/test_output.tsv"
    ARGS.validated = f"{top_dir}/6_out_hek293_glori_validated.pickle"
    ARGS.input_bed = f"{top_dir}/5_out_lifted/5_out_HEK293_SGNex_inosine_m6A_pileup_filter0.5_mod0.5.bed_cheui_like.bed_lifted"
    ARGS.output = f"{top_dir}/test_output.tsv"
    ARGS.aggregate = "avg"
    ARGS.metric = "rate"
    ARGS.discard = None
    ARGS.min_stoich = 0.8

    st = time.time()

    validated = ARGS.validated
    predicted_path = ARGS.input_bed
    out_path = ARGS.output

    if ARGS.metric not in ["model2","rate"]:
        raise ("--metric must be either model2 or rate")

    if ARGS.aggregate not in ["max","avg"]:
        raise ("--aggregate must be either max or avg")

    if ARGS.discard:
        with open(ARGS.discard,"rb") as p:
            discard_set = pickle.load(p)
    else:
        discard_set=set()

    with open(validated,"rb") as p:
        validated_set = pickle.load(p)

    p_lst_val, p_lst_all = [], [] # TODO: Not necessary
    validated_tested = 0 # TODO: Not necessary

    # Parse site level info to and store site thesholds
    preds_dct_all, preds_dct_val = parse_sites(predicted_path, validated_set, ARGS.metric, ARGS.min_stoich, ARGS.aggregate, discard_set)

    # Get list of probabilities
    p_lst_all = get_probs_list(preds_dct_all, ARGS.aggregate)
    p_lst_val = get_probs_list(preds_dct_val, ARGS.aggregate)

    # Initialise thresholds
    thresholds = define_thresholds()

    # Two pointers to get n predicted at thresholds for validated sites
    out_vals_all = get_n_predicted(thresholds, p_lst_all)
    out_vals_val = get_n_predicted(thresholds, p_lst_val)

    # Write output
    write_output(out_path, out_vals_val, out_vals_all)

    print("all done in ", time.time() - st, "seconds") # TODO: Convert to f-string

if __name__ == "__main__":
    # main()
    cProfile.run("main()", sort = "cumtime")
