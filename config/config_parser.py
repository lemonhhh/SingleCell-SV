import configargparse


def config_parser():
    parser = configargparse.ArgumentParser()
    # ---------- Global Setting Parameters --------- #
    # necessary for configargparse
    parser.add_argument('--config', is_config_file=True, help='config file path')
    # str example
    parser.add_argument("--expname", type=str, help='experiment name')
    parser.add_argument("--basedir", type=str, default='./logs/', help='where to store ckpts and logs')
    parser.add_argument("--datadir", type=str, default='.', help='input data directory')
    parser.add_argument("--labeldir", type=str, default='.', help='input label directory')

    # --------- training options --------- #
    # int example
    parser.add_argument("--netdepth", type=int, default=8, help='layers in coarse network')
    parser.add_argument("--netwidth", type=int, default=256, help='channels per layer in coarse network')
    # float example
    parser.add_argument("--lrate", type=float, default=5e-4, help='learning rate')
    parser.add_argument("--lrate_decay", type=float, default=0.1,
                        help='decay learning rate by a factor every specified number of steps')
    # batch size
    parser.add_argument("--N_rand", type=int, default=32 * 32 * 2,
                        help='batch size (number of random rays per gradient step)')
    parser.add_argument("--chunk", type=int, default=1024 * 8,
                        help='number of rays processed in parallel, decrease if running out of memory')
    # iterations
    parser.add_argument("--N_iters", type=int, default=250001, help='number of iterations')
    parser.add_argument("--no_reload", action='store_true',
                        help='do not reload weights from saved ckpt')
    parser.add_argument("--save_freq", type=int, default=5000, help='epoch frequency to save checkpoints')

    # ---------- rendering options ---------- #
    # store example, store_true denotes True if use_viewdirs exsists in the command line else False
    parser.add_argument("--use_viewdirs", action='store_true', help='use full 5D input instead of 3D')
    parser.add_argument("--N_samples", type=int, default=64, help='number of coarse samples per ray')
    parser.add_argument("--N_importance", type=int, default=0, help='number of additional fine samples per ray')
    parser.add_argument("--perturb", type=float, default=1.,
                        help='set to 0. for no jitter, 1. for jitter. Jitter for sampled t')
    parser.add_argument("--raw_noise_std", type=float, default=0.,
                        help='std dev of noise added to regularize sigma_a output, 1e0 recommended. Avoid floater')
    parser.add_argument("--multires", type=int, default=10,
                        help='log2 of max freq for positional encoding (3D location)')
    parser.add_argument("--multires_views", type=int, default=4,
                        help='log2 of max freq for positional encoding (2D direction)')
    parser.add_argument("--render_only", action='store_true',
                        help='do not optimize, reload weights and render out render_poses path')
    parser.add_argument("--render_test", action='store_true',
                        help='render the test set instead of render_poses path')

    # ----------- dataset options ------------ #
    parser.add_argument("--dataset_type", type=str, default='blender',
                        help='options: llff / blender / deepvoxels')
    # blender flags
    parser.add_argument("--white_bkgd", action='store_true',  # white background color for blender since alpha = 0
                        help='set to render synthetic data on a white bkgd (always use for dvoxels)')
    parser.add_argument("--half_res", action='store_true',
                        help='load blender synthetic data at 400x400 instead of 800x800')

    return parser


if __name__ == "__main__":
    parser = config_parser()
    args = parser.parse_args()

    print("debug")
