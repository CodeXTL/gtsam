import numpy as np
import gtsam
from gtsam import Values, NonlinearFactorGraph, LevenbergMarquardtOptimizer, noiseModel

# Import your custom SL4 class and bindings (assumed already wrapped)
from gtsam import SL4, PriorFactorSL4, BetweenFactorSL4


def main():
    # Random noise generator
    rng = np.random.default_rng(seed=42)

    def random_noise_vector(dim=15):
        return rng.uniform(low=-0.1, high=0.1, size=dim)

    # Create the factor graph
    graph = NonlinearFactorGraph()

    # Prior factor at node 1
    prior_sigmas = np.full(15, 0.3)
    prior_noise = noiseModel.Diagonal.Sigmas(prior_sigmas)
    graph.add(PriorFactorSL4(1, SL4(np.eye(4)), prior_noise))

    # # Between factor noise
    model_sigmas = np.full(15, 0.2)
    model = noiseModel.Diagonal.Sigmas(model_sigmas)

    # Define SL4 transformations (homography matrices)
    H12 = np.array(
        [
            [1.0, 0.1, 0.0, 2.0],
            [0.0, 1.0, 0.0, 3.0],
            [0.0, 0.0, 1.0, 5.0],
            [0.001, 0.002, 0.0, 1.0],
        ]
    )

    H23 = np.array(
        [
            [0.9, 0.2, 0.0, 1.5],
            [0.1, 1.1, 0.0, -2.0],
            [0.0, 0.0, 0.8, 4.0],
            [0.002, 0.003, 0.0005, 1.0],
        ]
    )

    H34 = np.array(
        [
            [1.05, -0.1, 0.0, 3.0],
            [0.2, 0.95, 0.0, 1.0],
            [0.0, 0.0, 0.9, 2.5],
            [0.0015, -0.001, 0.0003, 1.0],
        ]
    )

    H45 = np.array(
        [
            [0.98, 0.05, 0.0, -1.0],
            [-0.05, 1.02, 0.0, 2.0],
            [0.0, 0.0, 1.1, 0.5],
            [0.0008, 0.0015, -0.0002, 1.0],
        ]
    )

    H52 = np.linalg.inv(H23 @ H34 @ H45)

    # Create SL4 objects
    H12_SL4 = SL4(H12)
    H23_SL4 = SL4(H23)
    H34_SL4 = SL4(H34)
    H45_SL4 = SL4(H45)
    H52_SL4 = SL4(H52)
    H1 = SL4(np.eye(4))

    # Ground-truth poses
    H2 = H1.compose(H12_SL4)
    H3 = H2.compose(H23_SL4)
    H4 = H3.compose(H34_SL4)
    H5 = H4.compose(H45_SL4)

    gt_poses = [H1, H2, H3, H4, H5]

    # Double check the funtionalities
    # The below operation should be I_4x4
    print((H2.inverse() * H3) * H23_SL4.inverse())
    print("-------------------")
    # Each element should print the same thing
    print((H2.inverse() * H3))
    print(H23_SL4)
    print(H2.between(H3))
    print("-------------------")

    # Add odometry factors
    graph.add(BetweenFactorSL4(1, 2, H12_SL4, model))
    graph.add(BetweenFactorSL4(2, 3, H23_SL4, model))
    graph.add(BetweenFactorSL4(3, 4, H34_SL4, model))
    graph.add(BetweenFactorSL4(4, 5, H45_SL4, model))
    graph.add(BetweenFactorSL4(5, 2, H52_SL4, model))  # loop closure

    # Create initial estimate by perturbing GT
    initial = Values()
    for i, gt_pose in enumerate(gt_poses, 1):
        noise = random_noise_vector()
        noisy_pose = gt_pose.compose(SL4.Expmap(noise))
        initial.insert(i, noisy_pose)

    # Optimize
    params = gtsam.LevenbergMarquardtParams()
    params.setRelativeErrorTol(1e-5)
    params.setMaxIterations(100)
    optimizer = LevenbergMarquardtOptimizer(graph, initial, params)
    result = optimizer.optimize()

    print("Final result:")
    result.print()

    # Check accuracy
    for i, gt_pose in enumerate(gt_poses, 1):
        opt_pose = result.atSL4(i)
        if not gt_pose.equals(opt_pose, 1e-5):
            print(f"\033[1;31mPose {i} is outside tolerance!")
    print("\033[1;32mSuccessfully optimized!\033[0m")


if __name__ == "__main__":
    main()
