import numpy as np
# import cv2 # OpenCV is used for a direct comparison
import matplotlib.transforms as mtransforms

def get_affine_transform_numpy(src_p1, src_p2, src_p3, dst_p1, dst_p2, dst_p3):
    """
    Calculates the 2x3 affine transformation matrix from 3 pairs of
    corresponding points using NumPy.

    Args:
        src_p1 (list/tuple): Coordinates of the first source point [x, y].
        src_p2 (list/tuple): Coordinates of the second source point [x, y].
        src_p3 (list/tuple): Coordinates of the third source point [x, y].
        dst_p1 (list/tuple): Coordinates of the first destination point [x, y].
        dst_p2 (list/tuple): Coordinates of the second destination point [x, y].
        dst_p3 (list/tuple): Coordinates of the third destination point [x, y].

    Returns:
        np.ndarray: The 2x3 affine transformation matrix, or None if the
                    source points are collinear.
    """
    # Assemble the points into numpy arrays
    src_points = np.array([src_p1, src_p2, src_p3], dtype=np.float32)
    dst_points = np.array([dst_p1, dst_p2, dst_p3], dtype=np.float32)

    # We need to solve a system of linear equations.
    # The transformation is:
    # x' = a*x + b*y + c
    # y' = d*x + e*y + f
    #
    # For 3 points, we can set up a matrix equation:
    # [x1, y1, 1] [a]   [x1']
    # [x2, y2, 1] [b] = [x2']
    # [x3, y3, 1] [c]   [x3']
    #
    # And similarly for d, e, f:
    # [x1, y1, 1] [d]   [y1']
    # [x2, y2, 1] [e] = [y2']
    # [x3, y3, 1] [f]   [y3']

    # 1. Construct the coefficient matrix from the source points.
    # We add a column of ones to account for the translation (c and f).
    src_matrix = np.hstack([src_points, np.ones((3, 1))])

    # 2. Check for collinearity. If the determinant is zero, the points
    # are on a line and we can't find a unique solution.
    if np.linalg.det(src_matrix) == 0:
        print("Error: Source points are collinear. Cannot compute a unique affine transform.")
        return None

    # 3. Separate the destination points into x and y coordinates.
    dst_x = dst_points[:, 0]
    dst_y = dst_points[:, 1]

    # 4. Solve the system of linear equations for the two rows of the
    # affine matrix. np.linalg.solve is a robust way to do this.
    try:
        # Solve for the first row of the transform matrix (a, b, c)
        transform_row1 = np.linalg.solve(src_matrix, dst_x)
        # Solve for the second row of the transform matrix (d, e, f)
        transform_row2 = np.linalg.solve(src_matrix, dst_y)
    except np.linalg.LinAlgError:
        print("Error: Could not solve the system of equations. Check if points are valid.")
        return None

    # 5. Combine the two rows to form the 2x3 affine transformation matrix.
    transform_matrix = np.vstack([transform_row1, transform_row2])

    return transform_matrix

def get_affine_transform_opencv(src_p1, src_p2, src_p3, dst_p1, dst_p2, dst_p3):
    """
    Calculates the 2x3 affine transformation matrix using OpenCV's
    built-in function.

    Args:
        src_p1 (list/tuple): Coordinates of the first source point [x, y].
        src_p2 (list/tuple): Coordinates of the second source point [x, y].
        src_p3 (list/tuple): Coordinates of the third source point [x, y].
        dst_p1 (list/tuple): Coordinates of the first destination point [x, y].
        dst_p2 (list/tuple): Coordinates of the second destination point [x, y].
        dst_p3 (list/tuple): Coordinates of the third destination point [x, y].

    Returns:
        np.ndarray: The 2x3 affine transformation matrix.
    """
    # OpenCV requires the points to be of type float32.
    src_float = np.float32([src_p1, src_p2, src_p3])
    dst_float = np.float32([dst_p1, dst_p2, dst_p3])

    # cv2.getAffineTransform handles all the calculations internally.
    return cv2.getAffineTransform(src_float, dst_float)


def numpy_to_matplotlib_transform(numpy_matrix):
    """
    Converts a 2x3 NumPy/OpenCV affine matrix to a
    matplotlib.transforms.Affine2D object.

    Args:
        numpy_matrix (np.ndarray): A 2x3 affine transformation matrix.

    Returns:
        matplotlib.transforms.Affine2D: The corresponding Matplotlib transform.
    """
    # Matplotlib's Affine2D requires a 3x3 matrix.
    # We convert our 2x3 matrix to a 3x3 matrix by adding the [0, 0, 1] row.
    mpl_matrix = np.vstack([numpy_matrix, [0, 0, 1]])
    return mtransforms.Affine2D(mpl_matrix)


if __name__ == '__main__':
    # --- Define the 3 pairs of corresponding points ---
    # Points in the source plane
    p1_src = [50, 50]
    p2_src = [200, 50]
    p3_src = [50, 200]

    # Corresponding points in the destination plane
    # Let's say the transformation involves rotation and translation
    p1_dst = [100, 100]
    p2_dst = [236.6, 150.0]
    p3_dst = [50, 236.6]


    print("--- Calculating Affine Transform using NumPy ---")
    M_numpy = get_affine_transform_numpy(p1_src, p2_src, p3_src, p1_dst, p2_dst, p3_dst)
    if M_numpy is not None:
        print("Transformation Matrix:\n", M_numpy)

    print("\n" + "="*50 + "\n")

    print("--- Calculating Affine Transform using OpenCV ---")
    M_opencv = get_affine_transform_opencv(p1_src, p2_src, p3_src, p1_dst, p2_dst, p3_dst)
    print("Transformation Matrix:\n", M_opencv)

    # --- Verification ---
    if M_numpy is not None:
        print("\n" + "="*50 + "\n")
        print("--- Verifying the transformation ---")
        # Let's transform one of the source points and see if it matches the destination
        test_point_src = np.array([p1_src[0], p1_src[1], 1]) # Use homogeneous coordinates [x, y, 1]
        
        # Apply the transformation matrix
        # The result is [x', y']
        transformed_point = M_numpy @ test_point_src
        
        print(f"Transforming source point {p1_src}...")
        print(f"Calculated destination: {transformed_point}")
        print(f"Expected destination:   {p1_dst}")
        
        # Check if the calculated and expected points are close
        assert np.allclose(transformed_point, p1_dst)
        print("\nVerification successful!")
