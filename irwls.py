# Python Class for Iterative Weighted Least Squares in LDScore Regression
# Usage (single file): python3 irwls.py -d compare_em_ldsc/test/summary_data_h0.7_s1_r1_03151641235877.csv -o compare_em_ldsc/test -N 2000
# Usage (folder): python3 irwls.py -d compare_em_ldsc/test -o compare_em_ldsc/test -N 2000
# Save the log file in output_path/irwls.log and the result in output_path/irwls_results.csv
# -pd can be a fold or a file, if it is a fold, the program will run all files that contain "summary_data" in the fold

import numpy as np
import pandas as pd
import logging, time, argparse, os


class IRLS:
    """
    Python Class for Iterative Weighted Least Squares in LDScore Regression
    Parameters:
    :param ldscore: numpy array of shape (p, 1)
    :param zsq: squre of z-score: numpy array of shape (p, 1)
    :param weights: numpy array of shape (p, )
    Usage:
    ```
    from irwls import IRLS
    irwls = IRLS(ldscore, zsq, N)
    irwls.regression()
    reg_intercept = irwls.get_intercept()
    reg_h2 = irwls.get_coefficients()
    ```
    """

    def __init__(self, ldscore, zsq, N, weights=None):
        # initialize the class
        self.ldscore = ldscore
        self.zsq = zsq
        if weights is not None:
            self.weights = weights
        else:
            self.weights = np.ones_like(zsq)
        self.p = zsq.shape[0]
        self.N = N

    def _cal_initial_coef(self, max_int=1):
        """
        Estimate the initial regression coefficient with mean
        """
        zsq_sort = np.sort(self.zsq, axis=0)  # default: ascending
        intercept = np.mean(zsq_sort[: int(self.p * 0.95)])
        self.intercept = np.min([intercept, max_int])

        ldsc_mean = np.mean(self.ldscore, axis=0)
        zsq_mean = np.mean(self.zsq, axis=0)
        # self.beta = (zsq_mean - intercept) / ldsc_mean / self.N
        self.beta = (zsq_mean - self.intercept) / ldsc_mean / self.N * self.p
        # print("Initial intercept = %.4f, beta = %.8f" % (intercept, self.beta))

    def predict(self, beta=None, intercept=None):
        """
        Return the predicted value with given regression coefficient, or the current one in the class
        """
        if beta is None and intercept is None:
            if self.beta is not None and self.intercept is not None:
                beta, intercept = self.beta, self.intercept
            else:
                raise ValueError("No regression coefficient")
        return beta * self.ldscore * self.N / self.p + intercept

    def _update_weights(self, pred):
        """
        Update the weights with predicted value
        ref: ldsc/ldscore/regression.py-line: 498
        In ldsc.r: weight = 1 / (2 * pred ** 2)
        """
        self.weights /= 2 * pred**2  # /= var

    def regression(self, constraint_intercept=False, subtract=1):
        """
        Iterative Weighted Least Squares
        """
        self._cal_initial_coef()  # estimate raw coeffs with mean
        pred_raw = self.predict()
        self._update_weights(pred_raw)
        # get_coef in ldsc.r
        if constraint_intercept:
            zsq = self.zsq - subtract
            intercept = subtract
            ldscore = self.ldscore
        else:  # add new column of 1 to ldscore so that we can estimate the intercept
            ldscore = np.concatenate((np.ones((self.p, 1)), self.ldscore), axis=1)
            zsq = self.zsq
        w_ldsc = ldscore * np.sqrt(self.weights)
        w_zsq = zsq * np.sqrt(self.weights)
        # regression in ldsc.r
        # xtx = np.matmul(w_ldsc.T, w_ldsc)
        # xty = np.matmul(w_ldsc.T, w_zsq)
        # beta = np.linalg.solve(xtx, xty).flatten()
        beta = np.linalg.lstsq(w_ldsc, w_zsq, rcond=None)[0]
        if not constraint_intercept:
            self.intercept = beta[0].item()
            self.beta = beta[1].item() / self.N * self.p
        else:
            self.beta = beta.item() / self.N * self.p

    def get_intercept(self):
        """
        Return the intercept
        """
        return self.intercept

    def get_coefficients(self):
        """
        Return the coefficients
        """
        return self.beta


def get_parser():
    parser = argparse.ArgumentParser(description="Iterative Weighted Least Squares")
    parser.add_argument(
        "-d",
        "--data_path",
        type=str,
        default="compare_em_ldsc/test/summary_data_h0.6_s1_r1_03151640072308.csv",
        help="Path to the summary data",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
        default="compare_em_ldsc/test",
        help="Path to the output file",
    )
    parser.add_argument(
        "-N", "--N", type=int, default=1000, help="Number of individuals in LDSC"
    )
    return parser


if __name__ == "__main__":
    # set up parse arguments
    parser = get_parser()
    args = parser.parse_args()
    if not os.path.exists(args.output_path):
        # check if output path exists
        os.makedirs(args.output_path)
    # set up logging
    logger = logging.getLogger("irwls_logger")
    logger.setLevel(logging.INFO)
    log_path = args.output_path + "/log"
    if not os.path.exists(log_path):
        os.makedirs(log_path)
    fh = logging.FileHandler(log_path + "/irwls.log", mode="a")
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info("Start time: %s" % time.ctime())

    # get the data
    # if phenotype_data_path is a fold, get all files that contain "phenotype" in the fold
    summary_data_paths = []
    if os.path.isdir(args.data_path):
        summary_fold = args.data_path
        for root, dirs, files in os.walk(summary_fold):
            for file in files:
                if "summary_data" in file:
                    summary_data_paths.append(os.path.join(root, file))
    else:
        summary_data_paths.append(args.data_path)

    results = pd.DataFrame()

    # iterate through all files
    for summary_data_path in summary_data_paths:
        data = pd.read_csv(summary_data_path)
        ldscore = data["LDSCORE"].values.reshape(-1, 1)
        M = ldscore.shape[0]
        N = args.N
        zsq = data["Z"].values.reshape(-1, 1) ** 2

        # Get params from file name
        summary_params = (
            summary_data_path.split("summary_data_h")[1].split(".csv")[0].split("_")
        )
        hreal, sigma_beta, causal_rate = (
            float(summary_params[0][1:]),
            float(summary_params[1][1:]),
            float(summary_params[2][1:]),
        )
        created_time = summary_params[3]

        # fit the model
        irwls = IRLS(ldscore, zsq, N)
        irwls.regression()
        intercept = irwls.get_intercept()
        h_est_unfix = irwls.get_coefficients()
        irwls.regression(constraint_intercept=True, subtract=intercept)
        h_est_fix = irwls.get_coefficients()
        print("Intercept = %.4f, h = %.4f" % (intercept, h_est_unfix))
        logger.info("Intercept = %.4f, h = %.4f" % (intercept, h_est_unfix))
        print("h_fix = %.4f" % h_est_fix)

        result = pd.DataFrame(
            {
                "data": [created_time],
                "h": [hreal],
                "sigma_beta": [sigma_beta],
                "causal_rate": [causal_rate],
                "h_est_unfix": [round(h_est_unfix, 6)],
                "intercept": [round(intercept, 6)],
                "hest": [round(h_est_fix, 6)],
            }
        )
        results = pd.concat([results, result])

    # save the results to a file
    result_file = args.output_path + "/irwls_results.csv"
    with open(result_file, "a") as f:
        results.to_csv(f, header=f.tell() == 0, index=False)
    logger.info("End time: %s" % time.ctime())
