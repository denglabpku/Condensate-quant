import numpy as np
import os
import tifffile as tiff
import nd2
import random
import math
import torch
from torch.utils.data import Dataset
from skimage import io

def random_transform_one(input):
    """
    The function for data augmentation. Randomly select one method among five
    transformation methods (including rotation and flip) or do not use data
    augmentation.

    Args:
        input: the input patch before data augmentation
    Return:
        input: the input patch after data augmentation
    """
    p_trans = random.randrange(8)
    if p_trans == 0:  # no transformation
        input = input
    elif p_trans == 1:  # left rotate 90
        input = np.rot90(input, k=1, axes=(1, 2))
    elif p_trans == 2:  # left rotate 180
        input = np.rot90(input, k=2, axes=(1, 2))
    elif p_trans == 3:  # left rotate 270
        input = np.rot90(input, k=3, axes=(1, 2))
    elif p_trans == 4:  # horizontal flip
        input = input[:, :, ::-1]
    elif p_trans == 5:  # horizontal flip & left rotate 90
        input = input[:, :, ::-1]
        input = np.rot90(input, k=1, axes=(1, 2))
    elif p_trans == 6:  # horizontal flip & left rotate 180
        input = input[:, :, ::-1]
        input = np.rot90(input, k=2, axes=(1, 2))
    elif p_trans == 7:  # horizontal flip & left rotate 270
        input = input[:, :, ::-1]
        input = np.rot90(input, k=3, axes=(1, 2))
    return input

def random_transform(input, target):
    """
    The function for data augmentation. Randomly select one method among five
    transformation methods (including rotation and flip) or do not use data
    augmentation.

    Args:
        input, target : the input and target patch before data augmentation
    Return:
        input, target : the input and target patch after data augmentation
    """
    p_trans = random.randrange(8)
    if p_trans == 0:  # no transformation
        input = input
        target = target
    elif p_trans == 1:  # left rotate 90
        input = np.rot90(input, k=1, axes=(1, 2))
        target = np.rot90(target, k=1, axes=(1, 2))
    elif p_trans == 2:  # left rotate 180
        input = np.rot90(input, k=2, axes=(1, 2))
        target = np.rot90(target, k=2, axes=(1, 2))
    elif p_trans == 3:  # left rotate 270
        input = np.rot90(input, k=3, axes=(1, 2))
        target = np.rot90(target, k=3, axes=(1, 2))
    elif p_trans == 4:  # horizontal flip
        input = input[:, :, ::-1]
        target = target[:, :, ::-1]
    elif p_trans == 5:  # horizontal flip & left rotate 90
        input = input[:, :, ::-1]
        input = np.rot90(input, k=1, axes=(1, 2))
        target = target[:, :, ::-1]
        target = np.rot90(target, k=1, axes=(1, 2))
    elif p_trans == 6:  # horizontal flip & left rotate 180
        input = input[:, :, ::-1]
        input = np.rot90(input, k=2, axes=(1, 2))
        target = target[:, :, ::-1]
        target = np.rot90(target, k=2, axes=(1, 2))
    elif p_trans == 7:  # horizontal flip & left rotate 270
        input = input[:, :, ::-1]
        input = np.rot90(input, k=3, axes=(1, 2))
        target = target[:, :, ::-1]
        target = np.rot90(target, k=3, axes=(1, 2))
    return input, target

def n2v_mask_3d(volume, masking_ratio=0.05, neighbor_radius_xy=5, neighbor_radius_z=1):
    """
    3D Noise2Void blind-spot masking

    Args:
        volume : numpy array [D, H, W], the original input 3D patch (float or uint16, etc.)
        masking_ratio : float, the fraction of voxels to mask
        neighbor_radius_xy : int, the search radius for randomly selecting a neighbor voxel
        neighbor_radius_z : int, the search radius for randomly selecting a neighbor voxel

    Returns:
        input_masked : numpy array [D, H, W], with masked locations replaced by values from neighbors
        target       : numpy array [D, H, W], the original values (used as ground truth during training)
        mask         : numpy bool array [D, H, W], True indicates that the voxel is masked
    """
    D, H, W = volume.shape
    num_voxels = D * H * W
    k = max(1, int(masking_ratio * num_voxels))

    mask = np.zeros((D, H, W), dtype=bool)

    zs = np.random.randint(0, D, size=k)
    ys = np.random.randint(0, H, size=k)
    xs = np.random.randint(0, W, size=k)
    mask[zs, ys, xs] = True

    target = volume.copy()

    input_masked = volume.copy()

    for i in range(k):
        z, y, x = zs[i], ys[i], xs[i]

        dz = random.randint(-neighbor_radius_z, neighbor_radius_z)
        dy = random.randint(-neighbor_radius_xy, neighbor_radius_xy)
        dx = random.randint(-neighbor_radius_xy, neighbor_radius_xy)

        nz = np.clip(z + dz, 0, D - 1)
        ny = np.clip(y + dy, 0, H - 1)
        nx = np.clip(x + dx, 0, W - 1)

        if nz == z and ny == y and nx == x:
            nx = min(W - 1, nx + 1)
            ny = min(H - 1, ny + 1)

        input_masked[z, y, x] = volume[nz, ny, nx]

    return input_masked, target, mask

class trainset(Dataset):
    """
    Train set generator for pytorch training

    """

    def __init__(self, name_list, coordinate_list, noise_img_all, stack_index, masking_ratio=0.05, neighbor_radius_xy=5, neighbor_radius_z=1):
        self.name_list = name_list
        self.coordinate_list = coordinate_list
        self.noise_img_all = noise_img_all
        self.stack_index = stack_index
        self.masking_ratio = masking_ratio
        self.neighbor_radius_xy = neighbor_radius_xy
        self.neighbor_radius_z = neighbor_radius_z

    def __getitem__(self, index):
        """
        For temporal stacks with a small lateral size or short recording period, sub-stacks can be
        randomly cropped from the original stack to augment the training set according to the record
        coordinate. Then, interlaced frames of each sub-stack are extracted to form two 3D tiles.
        One of them serves as the input and the other serves as the target for network training
        Args:
            index : the index of 3D patchs used for training
        Return:
            input, target : the consecutive frames of the 3D noisy patch serve as the input and target of the network
        """
        stack_index = self.stack_index[index]
        noise_img = self.noise_img_all[stack_index]
        single_coordinate = self.coordinate_list[self.name_list[index]]
        init_h = single_coordinate['init_h']
        end_h = single_coordinate['end_h']
        init_w = single_coordinate['init_w']
        end_w = single_coordinate['end_w']
        init_s = single_coordinate['init_s']
        end_s = single_coordinate['end_s']
        input = noise_img[init_s:end_s, init_h:end_h, init_w:end_w]
        input = random_transform_one(input)

        input, target, mask = n2v_mask_3d(input, self.masking_ratio, self.neighbor_radius_xy, self.neighbor_radius_z)

        input = torch.from_numpy(np.expand_dims(input, 0).copy())
        target = torch.from_numpy(np.expand_dims(target, 0).copy())
        mask = torch.from_numpy(np.expand_dims(mask, 0).copy())
        return input, target, mask

    def __len__(self):
        return len(self.name_list)


class testset(Dataset):
    """
    Test set generator for pytorch inference

    """

    def __init__(self, name_list, coordinate_list, noise_img):
        self.name_list = name_list
        self.coordinate_list = coordinate_list
        self.noise_img = noise_img

    def __getitem__(self, index):
        """
        Generate the sub-stacks of the noisy image.
        Args:
            index : the index of 3D patch used for testing
        Return:
            noise_patch : the sub-stacks of the noisy image
            single_coordinate : the specific coordinate of sub-stacks in the noisy image for stitching all sub-stacks
        """
        single_coordinate = self.coordinate_list[self.name_list[index]]
        init_h = single_coordinate['init_h']
        end_h = single_coordinate['end_h']
        init_w = single_coordinate['init_w']
        end_w = single_coordinate['end_w']
        init_s = single_coordinate['init_s']
        end_s = single_coordinate['end_s']
        noise_patch = self.noise_img[init_s:end_s, init_h:end_h, init_w:end_w]
        noise_patch = torch.from_numpy(np.expand_dims(noise_patch, 0))
        return noise_patch, single_coordinate

    def __len__(self):
        return len(self.name_list)


def singlebatch_test_save(single_coordinate, output_image, raw_image):
    """
    Subtract overlapping regions (both the lateral and temporal overlaps) from the output sub-stacks (if the batch size equal to 1).

    Args:
        single_coordinate : the coordinate dict of the image
        output_image : the output sub-stack of the network
        raw_image : the noisy sub-stack
    Returns:
        output_patch : the output patch after subtract the overlapping regions
        raw_patch :  the raw patch after subtract the overlapping regions
        stack_start_ : the start coordinate of the patch in whole stack
        stack_end_ : the end coordinate of the patch in whole stack
    """
    stack_start_w = int(single_coordinate['stack_start_w'])
    stack_end_w = int(single_coordinate['stack_end_w'])
    patch_start_w = int(single_coordinate['patch_start_w'])
    patch_end_w = int(single_coordinate['patch_end_w'])

    stack_start_h = int(single_coordinate['stack_start_h'])
    stack_end_h = int(single_coordinate['stack_end_h'])
    patch_start_h = int(single_coordinate['patch_start_h'])
    patch_end_h = int(single_coordinate['patch_end_h'])

    stack_start_s = int(single_coordinate['stack_start_s'])
    stack_end_s = int(single_coordinate['stack_end_s'])
    patch_start_s = int(single_coordinate['patch_start_s'])
    patch_end_s = int(single_coordinate['patch_end_s'])

    output_patch = output_image[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
    raw_patch = raw_image[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
    return output_patch, raw_patch, stack_start_w, stack_end_w, stack_start_h, stack_end_h, stack_start_s, stack_end_s


def multibatch_test_save(single_coordinate, id, output_image, raw_image):
    """
    Subtract overlapping regions (both the lateral and temporal overlaps) from the output sub-stacks. (if the batch size larger than 1).

    Args:
        single_coordinate : the coordinate dict of the image
        output_image : the output sub-stack of the network
        raw_image : the noisy sub-stack
    Returns:
        output_patch : the output patch after subtract the overlapping regions
        raw_patch :  the raw patch after subtract the overlapping regions
        stack_start_ : the start coordinate of the patch in whole stack
        stack_end_ : the end coordinate of the patch in whole stack
    """
    stack_start_w_id = single_coordinate['stack_start_w'].numpy()
    stack_start_w = int(stack_start_w_id[id])
    stack_end_w_id = single_coordinate['stack_end_w'].numpy()
    stack_end_w = int(stack_end_w_id[id])
    patch_start_w_id = single_coordinate['patch_start_w'].numpy()
    patch_start_w = int(patch_start_w_id[id])
    patch_end_w_id = single_coordinate['patch_end_w'].numpy()
    patch_end_w = int(patch_end_w_id[id])

    stack_start_h_id = single_coordinate['stack_start_h'].numpy()
    stack_start_h = int(stack_start_h_id[id])
    stack_end_h_id = single_coordinate['stack_end_h'].numpy()
    stack_end_h = int(stack_end_h_id[id])
    patch_start_h_id = single_coordinate['patch_start_h'].numpy()
    patch_start_h = int(patch_start_h_id[id])
    patch_end_h_id = single_coordinate['patch_end_h'].numpy()
    patch_end_h = int(patch_end_h_id[id])

    stack_start_s_id = single_coordinate['stack_start_s'].numpy()
    stack_start_s = int(stack_start_s_id[id])
    stack_end_s_id = single_coordinate['stack_end_s'].numpy()
    stack_end_s = int(stack_end_s_id[id])
    patch_start_s_id = single_coordinate['patch_start_s'].numpy()
    patch_start_s = int(patch_start_s_id[id])
    patch_end_s_id = single_coordinate['patch_end_s'].numpy()
    patch_end_s = int(patch_end_s_id[id])

    output_image_id = output_image[id]
    raw_image_id = raw_image[id]
    output_patch = output_image_id[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]
    raw_patch = raw_image_id[patch_start_s:patch_end_s, patch_start_h:patch_end_h, patch_start_w:patch_end_w]

    return output_patch, raw_patch, stack_start_w, stack_end_w, stack_start_h, stack_end_h, stack_start_s, stack_end_s


def test_preprocess_lessMemoryNoTail_chooseOne(args, N):
    patch_y = args.patch_y
    patch_x = args.patch_x
    patch_t2 = args.patch_t
    gap_y = args.gap_y
    gap_x = args.gap_x
    gap_t2 = args.gap_t
    cut_w = (patch_x - gap_x) / 2
    cut_h = (patch_y - gap_y) / 2
    cut_s = (patch_t2 - gap_t2) / 2
    im_folder = args.datasets_path + '//' + args.datasets_folder

    name_list = []
    # train_raw = []
    coordinate_list = {}
    img_list = list(os.walk(im_folder, topdown=False))[-1][-1]
    img_list.sort()
    # print(img_list)

    im_name = img_list[N]

    im_dir = im_folder + '//' + im_name
    file_type = im_name[-3:]
    if file_type == 'tif':
        noise_im = tiff.imread(im_dir)
    elif file_type == 'nd2':
        noise_im = nd2.imread(im_dir)

    noise_im = tiff.imread(im_dir)
    # print('noise_im shape -----> ',noise_im.shape)
    # print('noise_im max -----> ',noise_im.max())
    # print('noise_im min -----> ',noise_im.min())
    if noise_im.shape[0] > args.test_datasize:
        noise_im = noise_im[0:args.test_datasize, :, :]
    noise_im = noise_im.astype(np.float32) / args.scale_factor
    # noise_im = (noise_im-noise_im.min()).astype(np.float32)/args.scale_factor

    whole_x = noise_im.shape[2]
    whole_y = noise_im.shape[1]
    whole_t = noise_im.shape[0]

    num_w = math.ceil((whole_x - patch_x + gap_x) / gap_x)
    num_h = math.ceil((whole_y - patch_y + gap_y) / gap_y)
    num_s = math.ceil((whole_t - patch_t2 + gap_t2) / gap_t2)
    # print('int((whole_y-patch_y+gap_y)/gap_y) -----> ',int((whole_y-patch_y+gap_y)/gap_y))
    # print('int((whole_x-patch_x+gap_x)/gap_x) -----> ',int((whole_x-patch_x+gap_x)/gap_x))
    # print('int((whole_t-patch_t2+gap_t2)/gap_t2) -----> ',int((whole_t-patch_t2+gap_t2)/gap_t2))
    for x in range(0, num_h):
        for y in range(0, num_w):
            for z in range(0, num_s):
                single_coordinate = {'init_h': 0, 'end_h': 0, 'init_w': 0, 'end_w': 0, 'init_s': 0, 'end_s': 0}
                if x != (num_h - 1):
                    init_h = gap_y * x
                    end_h = gap_y * x + patch_y
                elif x == (num_h - 1):
                    init_h = whole_y - patch_y
                    end_h = whole_y

                if y != (num_w - 1):
                    init_w = gap_x * y
                    end_w = gap_x * y + patch_x
                elif y == (num_w - 1):
                    init_w = whole_x - patch_x
                    end_w = whole_x

                if z != (num_s - 1):
                    init_s = gap_t2 * z
                    end_s = gap_t2 * z + patch_t2
                elif z == (num_s - 1):
                    init_s = whole_t - patch_t2
                    end_s = whole_t
                single_coordinate['init_h'] = init_h
                single_coordinate['end_h'] = end_h
                single_coordinate['init_w'] = init_w
                single_coordinate['end_w'] = end_w
                single_coordinate['init_s'] = init_s
                single_coordinate['end_s'] = end_s

                if y == 0:
                    single_coordinate['stack_start_w'] = y * gap_x
                    single_coordinate['stack_end_w'] = y * gap_x + patch_x - cut_w
                    single_coordinate['patch_start_w'] = 0
                    single_coordinate['patch_end_w'] = patch_x - cut_w
                elif y == num_w - 1:
                    single_coordinate['stack_start_w'] = whole_x - patch_x + cut_w
                    single_coordinate['stack_end_w'] = whole_x
                    single_coordinate['patch_start_w'] = cut_w
                    single_coordinate['patch_end_w'] = patch_x
                else:
                    single_coordinate['stack_start_w'] = y * gap_x + cut_w
                    single_coordinate['stack_end_w'] = y * gap_x + patch_x - cut_w
                    single_coordinate['patch_start_w'] = cut_w
                    single_coordinate['patch_end_w'] = patch_x - cut_w

                if x == 0:
                    single_coordinate['stack_start_h'] = x * gap_y
                    single_coordinate['stack_end_h'] = x * gap_y + patch_y - cut_h
                    single_coordinate['patch_start_h'] = 0
                    single_coordinate['patch_end_h'] = patch_y - cut_h
                elif x == num_h - 1:
                    single_coordinate['stack_start_h'] = whole_y - patch_y + cut_h
                    single_coordinate['stack_end_h'] = whole_y
                    single_coordinate['patch_start_h'] = cut_h
                    single_coordinate['patch_end_h'] = patch_y
                else:
                    single_coordinate['stack_start_h'] = x * gap_y + cut_h
                    single_coordinate['stack_end_h'] = x * gap_y + patch_y - cut_h
                    single_coordinate['patch_start_h'] = cut_h
                    single_coordinate['patch_end_h'] = patch_y - cut_h

                if z == 0:
                    single_coordinate['stack_start_s'] = z * gap_t2
                    single_coordinate['stack_end_s'] = z * gap_t2 + patch_t2 - cut_s
                    single_coordinate['patch_start_s'] = 0
                    single_coordinate['patch_end_s'] = patch_t2 - cut_s
                elif z == num_s - 1:
                    single_coordinate['stack_start_s'] = whole_t - patch_t2 + cut_s
                    single_coordinate['stack_end_s'] = whole_t
                    single_coordinate['patch_start_s'] = cut_s
                    single_coordinate['patch_end_s'] = patch_t2
                else:
                    single_coordinate['stack_start_s'] = z * gap_t2 + cut_s
                    single_coordinate['stack_end_s'] = z * gap_t2 + patch_t2 - cut_s
                    single_coordinate['patch_start_s'] = cut_s
                    single_coordinate['patch_end_s'] = patch_t2 - cut_s

                # noise_patch1 = noise_im[init_s:end_s,init_h:end_h,init_w:end_w]
                patch_name = args.datasets_folder + '_x' + str(x) + '_y' + str(y) + '_z' + str(z)
                # train_raw.append(noise_patch1.transpose(1,2,0))
                name_list.append(patch_name)
                # print(' single_coordinate -----> ',single_coordinate)
                coordinate_list[patch_name] = single_coordinate

    return name_list, noise_im, coordinate_list


def test_preprocess_chooseOne(args, img_id):
    """
    Choose one original noisy stack and partition it into thousands of 3D sub-stacks (patch) with the setting
    overlap factor in each dimension.

    Args:
        args : the train object containing input params for partition
        img_id : the id of the test image
    Returns:
        name_list : the coordinates of 3D patch are indexed by the patch name in name_list
        noise_im : the original noisy stacks
        coordinate_list : record the coordinate of 3D patch preparing for partition in whole stack
        im_name : the file name of the noisy stacks

    """

    patch_y = args.patch_y
    patch_x = args.patch_x
    patch_t = args.patch_t
    gap_y = args.gap_y
    gap_x = args.gap_x
    gap_t = 1 if args.patch_t == 1 else int(args.patch_t * (1 - args.overlap_factor))
    cut_w = (patch_x - gap_x) / 2
    cut_h = (patch_y - gap_y) / 2
    cut_s = (patch_t - gap_t) / 2
    im_folder = args.datasets_path

    name_list = []
    coordinate_list = {}
    img_list = list(os.walk(im_folder, topdown=False))[-1][-1]
    img_list.sort()

    im_name = img_list[img_id]


    im_dir = im_folder + '//' + im_name

    file_type = im_name[-3:]
    if file_type == 'tif':
        noise_im = tiff.imread(im_dir)
    elif file_type == 'nd2':
        noise_im = nd2.imread(im_dir)

    input_data_type = noise_im.dtype
    img_mean = noise_im.mean()
    # print('noise_im max -----> ',noise_im.max())
    # print('noise_im min -----> ',noise_im.min())
    if noise_im.shape[0] > args.test_datasize:
        noise_im = noise_im[0:args.test_datasize, :, :]
    if args.print_img_name:
       print('Testing image name -----> ', im_name)
       print('Testing image shape -----> ', noise_im.shape)
    # Minus mean before training
    noise_im = noise_im.astype(np.float32)/args.scale_factor
    noise_im = noise_im-img_mean
    # No preprocessing
    # noise_im = noise_im.astype(np.float32) / args.scale_factor
    # noise_im = (noise_im-noise_im.min()).astype(np.float32)/args.scale_factor

    whole_x = noise_im.shape[2]
    whole_y = noise_im.shape[1]
    whole_t = noise_im.shape[0]

    num_w = math.ceil((whole_x - patch_x + gap_x) / gap_x)
    num_h = math.ceil((whole_y - patch_y + gap_y) / gap_y)
    num_s = math.ceil((whole_t - patch_t + gap_t) / gap_t)
    # print('int((whole_y-patch_y+gap_y)/gap_y) -----> ',int((whole_y-patch_y+gap_y)/gap_y))
    # print('int((whole_x-patch_x+gap_x)/gap_x) -----> ',int((whole_x-patch_x+gap_x)/gap_x))
    # print('int((whole_t-patch_t2+gap_t2)/gap_t2) -----> ',int((whole_t-patch_t2+gap_t2)/gap_t2))
    for x in range(0, num_h):
        for y in range(0, num_w):
            for z in range(0, num_s):
                single_coordinate = {'init_h': 0, 'end_h': 0, 'init_w': 0, 'end_w': 0, 'init_s': 0, 'end_s': 0}
                if x != (num_h - 1):
                    init_h = gap_y * x
                    end_h = gap_y * x + patch_y
                elif x == (num_h - 1):
                    init_h = whole_y - patch_y
                    end_h = whole_y

                if y != (num_w - 1):
                    init_w = gap_x * y
                    end_w = gap_x * y + patch_x
                elif y == (num_w - 1):
                    init_w = whole_x - patch_x
                    end_w = whole_x

                if z != (num_s - 1):
                    init_s = gap_t * z
                    end_s = gap_t * z + patch_t
                elif z == (num_s - 1):
                    init_s = whole_t - patch_t
                    end_s = whole_t
                single_coordinate['init_h'] = init_h
                single_coordinate['end_h'] = end_h
                single_coordinate['init_w'] = init_w
                single_coordinate['end_w'] = end_w
                single_coordinate['init_s'] = init_s
                single_coordinate['end_s'] = end_s

                if y == 0:
                    single_coordinate['stack_start_w'] = y * gap_x
                    single_coordinate['stack_end_w'] = y * gap_x + patch_x - cut_w
                    single_coordinate['patch_start_w'] = 0
                    single_coordinate['patch_end_w'] = patch_x - cut_w
                elif y == num_w - 1:
                    single_coordinate['stack_start_w'] = whole_x - patch_x + cut_w
                    single_coordinate['stack_end_w'] = whole_x
                    single_coordinate['patch_start_w'] = cut_w
                    single_coordinate['patch_end_w'] = patch_x
                else:
                    single_coordinate['stack_start_w'] = y * gap_x + cut_w
                    single_coordinate['stack_end_w'] = y * gap_x + patch_x - cut_w
                    single_coordinate['patch_start_w'] = cut_w
                    single_coordinate['patch_end_w'] = patch_x - cut_w

                if x == 0:
                    single_coordinate['stack_start_h'] = x * gap_y
                    single_coordinate['stack_end_h'] = x * gap_y + patch_y - cut_h
                    single_coordinate['patch_start_h'] = 0
                    single_coordinate['patch_end_h'] = patch_y - cut_h
                elif x == num_h - 1:
                    single_coordinate['stack_start_h'] = whole_y - patch_y + cut_h
                    single_coordinate['stack_end_h'] = whole_y
                    single_coordinate['patch_start_h'] = cut_h
                    single_coordinate['patch_end_h'] = patch_y
                else:
                    single_coordinate['stack_start_h'] = x * gap_y + cut_h
                    single_coordinate['stack_end_h'] = x * gap_y + patch_y - cut_h
                    single_coordinate['patch_start_h'] = cut_h
                    single_coordinate['patch_end_h'] = patch_y - cut_h

                if z == 0:
                    single_coordinate['stack_start_s'] = z * gap_t
                    single_coordinate['stack_end_s'] = z * gap_t + patch_t - cut_s
                    single_coordinate['patch_start_s'] = 0
                    single_coordinate['patch_end_s'] = patch_t - cut_s
                elif z == num_s - 1:
                    single_coordinate['stack_start_s'] = whole_t - patch_t + cut_s
                    single_coordinate['stack_end_s'] = whole_t
                    single_coordinate['patch_start_s'] = cut_s
                    single_coordinate['patch_end_s'] = patch_t
                else:
                    single_coordinate['stack_start_s'] = z * gap_t + cut_s
                    single_coordinate['stack_end_s'] = z * gap_t + patch_t - cut_s
                    single_coordinate['patch_start_s'] = cut_s
                    single_coordinate['patch_end_s'] = patch_t - cut_s

                # noise_patch1 = noise_im[init_s:end_s,init_h:end_h,init_w:end_w]
                patch_name = args.datasets_name + '_x' + str(x) + '_y' + str(y) + '_z' + str(z)
                # train_raw.append(noise_patch1.transpose(1,2,0))
                name_list.append(patch_name)
                # print(' single_coordinate -----> ',single_coordinate)
                coordinate_list[patch_name] = single_coordinate

    return name_list, noise_im, coordinate_list, im_name, img_mean, input_data_type
