# Use EM to solve linear mixed model
# Usage: 
# python3 /home/wjiang49/UKBheight/lmm_em.py -pd /home/wjiang49/scratch/UKBsimdata/h0.6/phenotypeh0.6_s2_r0.01_0228213237.csv -gd /home/wjiang49/scratch/UKBsimdata/genotypes.csv -o /home/wjiang49/scratch/
# Save the log file in output_path/em.log and the result in output_path/em_results.csv

import numpy as np
import pandas as pd
import time, argparse, logging, os


def lmm_em(y, X, Z, tol=1e-3, max_iter=10, verbose=True):
    """
    Input:
    y: n x 1, response
    X: n x p, random effect
    Z: n x c, fixed effect
    tol: float, tolerance for convergence, default 1e-3
    max_iter: int, maximum number of iterations, default 10
    verbose: bool, print information, default True

    Output:
    likelihood_list: likelihood in each iteration
    omega_list: value x (iter + 1), record omega in each iteration
    sigma_beta2_list: record sigma_beta^2 in each iteration
    sigma_e2_list: record sigma_e^2 in each iteration
    beta_post_mean: posterior mean of beta

    Notation:
    omega: c x 1, fixed effect coefficient
    beta: p x 1, random effect coefficient

    Formula and prior distribution:
    y = Z * omega + X * beta + e
    e ~ N(0, sigma_e^2 * I_n)
    beta ~ N(0, sigma_beta^2 * I_p)

    Initialization:
    omega = (Z^T * Z)^-1 * Z^T * y
    sigma_beta^2 = var(y - Z * omega) / 2
    sigma_e^2 = var(y - Z * omega) / 2

    Implementation:
    Eignvalue decomposition:
        X^T * X = eigvecs_xtx * diag(eigvals_xtx) * eigvecs_xtx^T (Q: p x p)
        X * X^T = eigvecs_xxt * diag(eigvals_xxt) * eigvecs_xxt^T (tilde{Q}: n x n)

    1. For n >= p:
    E-step:
        Gamma = eigvecs_xtx * diag(1 / (eigvals_xtx / sigma_e^2 + \
            1 / sigma_beta^2)) * eigvecs_xtx^T
        mu = Gamma * X^T * (y - Z * omega) / sigma_e^2

    M-step:
        omega = (Z^T * Z)^-1 * Z^T * (y - X * mu)
        sigma_beta^2 = (sum sigma_e^2 * sigma_beta^2 / \
            (sigma_e^2 + sigma_beta^2 * eigvals_xtx)+ ||mu||^2) / p
        sigma_e^2 = (||y - Z * omega||^2 + sum(eigvals_xtx * sigma_beta^2 * sigma_e^2 / \
            (eigvals_xtx * sigma_beta^2 + sigma_e^2)) + \
            ||X * mu||^2 - 2 * (y - Z * omega)^T * X * mu) / n

    2. For n < p:
    E-step:
        Gamma = eigvecs_xxt * diag(1 / (eigvals_xxt / sigma_e^2 + \
            1 / sigma_beta^2)) * eigvecs_xxt^T
        mu = X^T * Gamma * (y - Z * omega) / sigma_e^2

    M-step:
        omega = (Z^T * Z)^-1 * Z^T * (y - X * mu)
        sigma_beta^2 = (sum sigma_e^2 * sigma_beta^2 / \
            (sigma_e^2 + sigma_beta^2 * eigvals_xxt) \
            + ||mu||^2 + (p - n) * sigma_beta^2) / p
        sigma_e^2 = (||y - Z * omega||^2 + sum(eigvals_xxt * sigma_beta^2 * sigma_e^2 / \
            (eigvals_xxt * sigma_beta^2 + sigma_e^2)) + \
            ||X * mu||^2 - 2 * (y - Z * omega)^T * X * mu) / n

    log-likelihood (n >= p and n < p):
        l = - n / 2 * log(2 * pi) - 1 / 2 * log(sum (sigma_beta^2 * eigvals_xxt + \
            sigma_e^2)) - (y - Z * omega)^T * eigvecs_xxt * \
            1 / diag(sigma_beta^2 * eigvals_xxt + sigma_e^2) * eigvecs_xxt^T * \
            (y - Z * omega) / 2
    """
    n, p = X.shape
    n_, c = Z.shape
    assert n == n_, "X and Z must have same number of rows"
    assert y.shape == (n, 1), "y must be a column vector with length n"

    # Choose the calculation method according to n and p
    n_geq_p = n >= p
    print("n = {}, p = {}, model: {}".format(n, p, "n >= p" if n_geq_p else "n < p"))

    # Initialize parameters
    mu = np.random.randn(p, 1)
    omega = np.linalg.inv(Z.T @ Z) @ Z.T @ y
    y_z_omega = y - Z @ omega
    sigma_beta2 = np.var(y_z_omega) / 2
    sigma_e2 = sigma_beta2

    # For accelerate calculation
    XXT = X @ X.T
    XTX = X.T @ X
    ZTZinvZT = np.linalg.inv(Z.T @ Z) @ Z.T

    # eigenvalue decomposition of XXT, needed in n < p & l
    eigvals_xxt, eigvecs_xxt = np.linalg.eigh(XXT)  # l, ~Q

    if n_geq_p:
        # eigenvalue decomposition of XTX, needed in n >= p
        eigvals_xtx, eigvecs_xtx = np.linalg.eigh(XTX)  # Q

        # n >= p
        def cal_mu():
            return (
                eigvecs_xtx
                @ (eigvecs_xtx.T @ (X.T @ y_z_omega / d_.reshape(-1, 1)))
                / sigma_e2
            )

        def cal_sigma_beta2():
            return (
                np.sum(sigma_e2 * sigma_beta2 / (sigma_e2 + sigma_beta2 * eigvals_xtx))
                + np.linalg.norm(mu) ** 2
            ) / p

        def cal_sigma_e2():
            return (
                np.linalg.norm(y_z_omega) ** 2
                + np.sum(
                    eigvals_xtx
                    * sigma_beta2
                    * sigma_e2
                    / (eigvals_xtx * sigma_beta2 + sigma_e2)
                )
                + np.linalg.norm(X @ mu) ** 2
                - 2 * (y_z_omega.T @ X @ mu)[0, 0]
            ) / n

    else:
        # n < p
        def cal_mu():
            return (
                X.T
                @ (eigvecs_xxt @ (eigvecs_xxt.T @ y_z_omega / d_.reshape(-1, 1)))
                / sigma_e2
            )

        def cal_sigma_beta2():
            return (
                np.sum(1 / d_) + np.linalg.norm(mu) ** 2 + (p - n) * sigma_beta2
            ) / p

        def cal_sigma_e2():
            return (
                np.linalg.norm(y_z_omega - X @ mu) ** 2 + np.sum(eigvals_xxt / d_)
            ) / n

    # log-likelihood
    likelihood_const = -n / 2 * np.log(2 * np.pi)

    def cal_likelihood():
        return (
            likelihood_const
            - np.sum(np.log(sigma_beta2 * eigvals_xxt + sigma_e2)) / 2
            - y_z_omega.T
            @ eigvecs_xxt
            @ np.diag(1 / (sigma_beta2 * eigvals_xxt + sigma_e2))
            @ eigvecs_xxt.T
            @ y_z_omega
            / 2
        )

    # Record parameters in each iteration
    max_iter += 1
    likelihood_list = np.zeros(max_iter)
    omega_list = np.zeros((c, max_iter))
    sigma_beta2_list = np.zeros(max_iter)
    sigma_e2_list = np.zeros(max_iter)

    # Initialize record
    iter = 0
    omega_list[:, iter] = omega.flatten()
    likelihood_list[iter] = np.squeeze(cal_likelihood())
    sigma_beta2_list[iter] = sigma_beta2
    sigma_e2_list[iter] = sigma_e2

    # EM algorithm
    print("EM algorithm starts")
    convergence = False
    for iter in range(1, max_iter):
        # E step
        d_ = (
            eigvals_xxt / sigma_e2 + 1 / sigma_beta2
            if not n_geq_p
            else eigvals_xtx / sigma_e2 + 1 / sigma_beta2
        )
        mu = cal_mu()

        # M step
        omega = ZTZinvZT @ (y - X @ mu)
        y_z_omega = y - Z @ omega
        sigma_beta2 = cal_sigma_beta2()
        sigma_e2 = cal_sigma_e2()

        # Calculate likelihood
        likelihood_list[iter] = np.squeeze(cal_likelihood())

        # Record parameters
        omega_list[:, iter] = omega.flatten()
        sigma_beta2_list[iter] = sigma_beta2
        sigma_e2_list[iter] = sigma_e2

        # Check convergence
        if np.abs(likelihood_list[iter] - likelihood_list[iter - 1]) < tol:
            convergence = True
            break

        # Print process information
        if verbose and iter % 10 == 0:
            print(
                "iter: {}, log-likelihood: {:.4e}".format(iter, likelihood_list[iter])
            )
            print("beta: {:.4e}".format(np.mean(mu)))
            print("sigma_beta^2: {:.4e}".format(sigma_beta2))
            print("sigma_e^2: {:.4e}".format(sigma_e2))
            print("--------------------------------------")

    # Algorithm summary
    beta_post_mean = np.mean(mu)
    resident = np.linalg.norm(y - Z @ omega - X @ mu) ** 2 / n
    if not convergence:
        print("EM algorithm does not converge within {} iterations".format(iter))
    else:
        print("EM algorithm converges after {} iterations".format(iter))
    print("sigma_beta^2 = {:.4e}".format(sigma_beta2))
    print("sigma_e^2 = {:.4e}".format(sigma_e2))
    print("beta_post_mean = {:.4e}".format(beta_post_mean))
    print("omega_mean = {:.4e}".format(np.mean(omega)))
    print("resident = {:.4e}".format(resident))

    return (
        likelihood_list[: iter + 1],
        omega_list[:, : iter + 1],
        sigma_beta2_list[: iter + 1],
        sigma_e2_list[: iter + 1],
        mu # beta_post
    )


def get_parser():
    parser = argparse.ArgumentParser(description="LMM EM Algorithm")
    parser.add_argument(
        "-pd",
        "--phenotype_data_path",
        default="phenotypeh0.8.csv",
        help="Path to the summary data"
    )
    parser.add_argument(
        "-gd",
        "--genotype_data_path",
        type=str,
        default="genotypes_partial.csv",
        help="Path to the summary data"
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
        default="em",
        help="Path to the output file"
    )
    return parser

if __name__ == "__main__":
    # set up logging and parse arguments
    parser = get_parser()
    args = parser.parse_args()
    if not os.path.exists(args.output_path):
        # check if output path exists
        os.makedirs(args.output_path)
    logger = logging.getLogger('em_logger')
    logger.setLevel(logging.INFO)
    # fh = logging.FileHandler(args.output_path+"/em.log", mode='a')
    # fh.setLevel(logging.INFO)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)
    logger.info("Start time: %s" % time.ctime())
    print("Start time: %s" % time.ctime())

    # load data
    genotypes_path = args.genotype_data_path
    genotypes = pd.read_csv(genotypes_path, header=0).values

    # get all files that contain "phenotype" in the fold
    phenotype_fold = args.phenotype_data_path
    phenotype_paths = []
    for root, dirs, files in os.walk(phenotype_fold):
        for file in files:
            if "phenotype" in file:
                phenotype_paths.append(os.path.join(root, file))

    for phenotype_path in phenotype_paths:
        phenotype = pd.read_csv(phenotype_path, header=0).values

        y = phenotype.reshape(-1, 1)
        X = genotypes
        Z = np.ones((y.shape[0], 1))

        n, p = X.shape

        # run EM algorithm
        start_time = time.time()
        # run EM algorithm
        (
            likelihood_list,
            omega_list,
            sigma_beta2_list,
            sigma_e2_list,
            beta_post
        ) = lmm_em(y, X, Z, tol=1e-6, max_iter=2000, verbose=False)
        end_time = time.time()
        print(
            "Run time: %d min %.2f s"
            % ((end_time - start_time) // 60, (end_time - start_time) % 60)
        )
        sigma_beta2 = sigma_beta2_list[-1]
        sigma_e2 = sigma_e2_list[-1]
        # logger.info("sigma_beta^2 = {:.4e}".format(sigma_beta2))
        # logger.info("sigma_e^2 = {:.4e}".format(sigma_e2))
        # logger.info("beta_post_mean = {:.4e}".format(np.mean(beta_post)))
        # logger.info("omega_mean = {:.4e}".format(np.mean(omega_list[:, -1])))

        # Get params from file name
        # file_name: phenotypeh${heritability}_s${sigma_beta}_r${causal_rate}_%m%d%H%M%S.csv
        phenotype_params = phenotype_path.split("phenotypeh")[1].split(".csv")[0].split("_")
        h, sigma_beta, causal_rate = float(phenotype_params[0][1:]), float(phenotype_params[1][1:]), float(phenotype_params[2][1:])
        created_time = phenotype_params[3]

        # Compare with real h
        var_xbeta = p * sigma_beta2
        h_calc = var_xbeta / (var_xbeta + sigma_e2)
        logger.info("Real h = %.4f, calculated h = %.4f" % (h, h_calc))
        print("Real h2g = %.4f, calculated h2g = %.4f" % (h, h_calc))

        results = pd.DataFrame({
            "data": [created_time],
            "h": [h],
            "sigma_beta": [sigma_beta],
            "causal_rate": [causal_rate],
            "hest": [round(h_calc, 6)],
            "sigma_beta_est": [round(np.sqrt(sigma_beta2), 6)]
        })

        # Save result
        result_file = args.output_path + "/em_results.csv"
        with open(result_file, "a") as f:
            results.to_csv(f, header=f.tell()==0, index=False)

        logger.info("End time: %s\n\n" % time.ctime())
