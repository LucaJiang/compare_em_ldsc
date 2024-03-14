# Python Class for Iterative Weighted Least Squares

import numpy as np
import pandas as pd
import logging, time, argparse
# Usage: python3 irwls.py -d /home/wjiang49/scratch/summarydata/summary_datah0.6.csv -o /home/wjiang49/scratch/summarydata/h2g


class IRLS:
    """
    Python Class for Iterative Weighted Least Squares
    Parameters:
    :param X: numpy array of shape (n_samples, 1)
    :param y: numpy array of shape (n_samples, 1)
    :param fix_intercept: bool, whether to fix the intercept to 0, default False
    :param weights: numpy array of shape (n_samples, )

    Usage:
    ```
    from irwls import IRLS
    irwls = IRLS(X, y)
    irwls.regression()
    reg_intercept = irwls.get_intercept()
    reg_coefficients = irwls.get_coefficients()
    ```
    """

    def __init__(
        self,
        X,
        y,
        fix_intercept=False,
        weights=None,
        max_iter=1,
    ):
        """
        :param X: numpy array of shape (n_samples, n_features)
        :param y: numpy array of shape (n_samples, )
        :param fix_intercept: bool, whether to fix the intercept to 0, default False
        :param weights: numpy array of shape (n_samples, )
        :param max_iter: int, maximum number of iterations
        :param tol: float, tolerance for stopping criteria
        """
        if not fix_intercept:
            self.X = np.concatenate([np.ones_like(X), X], axis=1)
        else:
            self.X = X
        self.y = y
        self.fix_intercept = fix_intercept
        if weights is not None:
            self.weights = weights
        else:
            self.weights = np.ones_like(y)
        self.max_iter = max_iter
        self.p = 1

    def _get_initial_beta(self):
        """
        Initialize the regression coefficient with one
        """
        # return np.ones_like(self.X[0])
        return np.linalg.pinv(self.X) @ self.y

    def _get_residuals(self, X, y, beta):
        return y - X @ beta

    def _get_weights_r(self, residuals):
        """
        Calculate the weight vector with residuals
        ref: initial weighted method of irwls
        """
        weights = np.abs(residuals) ** (self.p - 2)
        weights /= weights.sum()
        return weights

    def _get_weights_yhat(self, pred):
        """
        Calculate the weight vector with residuals
        ref: ldsc/ldscore/regression.py-line: 498
        In ldsc.r: weight = 1 / (2 * pred ** 2)
        """
        weights = pred ** (-1 / 2) / 2 / self.weights
        return weights

    def predict(self, X):
        return X @ self.beta

    def _get_beta(self, X, y, weights):
        """
        Update the regression coefficient with weighted least squares
        return: numpy array of shape (n_features, 1)
        """
        return np.linalg.lstsq(
            np.multiply(X, weights), np.multiply(y, weights), rcond=None
        )[0].reshape(-1, 1)

    def _check_convergence(self, residuals, new_residuals):
        return np.linalg.norm(residuals - new_residuals) < self.tol

    def regression(self):
        """
        Fit the model and return the regression coefficient
        """
        # Initialize the regression coefficient and residuals
        self.beta = self._get_initial_beta()
        # Iterate
        for _ in range(self.max_iter):
            # residuals = self._get_residuals(self.X, self.y, beta)
            # Calculate the weight vector
            weights = self._get_weights_r(self.predict(self.X))
            # Update the regression coefficient with weighted least squares
            self.beta = self._get_beta(self.X, self.y, weights)
        # Save results
        # self.residuals = residuals

    def get_intercept(self):
        """
        Return the intercept
        """
        return self.beta[0].item()

    def get_coefficients(self):
        """
        Return the coefficients
        """
        if self.fix_intercept:
            return self.beta.item()
        return self.beta[1].item()

def save_results(path, intercept, h2, h2_fix, real_h2g):
    if path is not None:
        with open(path, "a") as f:
            f.write("%.4f %.4f %.4f %.4f\n" % (intercept, h2, h2_fix, real_h2g))
    else:
        df = pd.DataFrame({"intercept": intercept, "h2": h2, "h2_fix": h2_fix, "real_h2g": real_h2g})
        df.to_csv(path, index=False, sep="\t", header=True)


def get_parser():
    parser = argparse.ArgumentParser(description="Iterative Weighted Least Squares")
    parser.add_argument(
        "-d",
        "--data_path",
        type=str,
        default="summary_datah0.6.csv",
        help="Path to the summary data"
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
        default="h2g",
        help="Path to the output file"
    )
    parser.add_argument(
        "-N",
        "--N",
        type=int,
        default=50000,
        help="Number of individuals in LDSC"
    )
    parser.add_argument(
        "-M",
        "--M",
        type=int,
        default=44717,
        help="Number of SNPs in LDSC"
    )

    return parser

#! TEST
if __name__ == "__main__":
    # set up logging and parse arguments
    parser = get_parser()
    args = parser.parse_args()
    logger = logging.getLogger('irwls_logger')
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(args.output_path+".log", mode='a')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.info("Start time: %s" % time.ctime())

    data = pd.read_csv(args.data_path)
    # get real h from name of the file
    real_h2g = float(args.data_path.split("summary_datah")[1].split(".csv")[0])
    logging.info("Using data: %s with real h_g^2 = %.4f" % (args.data_path, real_h2g))

    y = data["Z"].values.reshape(-1, 1)**2
    l2 = data["LDSCORE"].values.reshape(-1, 1)
    # x = l2 * 61220 / 853604  # l2 * N / M
    x = l2 * args.N / args.M  # l2 * N / M
    irwls = IRLS(x, y, weights=l2)
    irwls.regression()
    # print("Intercept = %.4f" % irwls.get_intercept())
    # print("h_g^2 = %.4f" % irwls.get_coefficients())
    intercept = irwls.get_intercept()
    h2g = irwls.get_coefficients()
    logger.info("Intercept = %.4f" % intercept)
    logger.info("h_g^2 = %.4f" % h2g)

    irwls_fix = IRLS(x, y - intercept, fix_intercept=True)
    irwls_fix.regression()
    # print("h_g^2 fix = %.4f" % irwls_fix.get_coefficients())
    h2g_fix = irwls_fix.get_coefficients()
    logger.info("h_g^2 fix = %.4f" % h2g_fix)
    save_results(args.output_path+".csv", intercept, h2g, h2g_fix, real_h2g)
    logger.info("End time: %s" % time.ctime())

