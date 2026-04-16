#include "identification.h"

list_node* base_table[N_STARS][CONFIG_SIZE - 1];


star_tracker_error list_insert(list_node **first, int val) {
    if (first == NULL) {
        return incorrect_input;
    }
    list_node *new = malloc(sizeof(list_node));
    if (new == NULL) {
        return memory_error;
    }
    new->next = *first;
    new->val = val;
    *first = new;
    return ok;
}


star_tracker_error list_clear(list_node **first){
    if (first == NULL) {
        return incorrect_input;
    }
    list_node *tmp;
    while (*first != NULL) {
        tmp = *first;
        *first = (*first)->next;
        free(tmp);
    }
    return ok;
}


star_tracker_error set_dist(int star1, int star2, int N_stars, double *dist_mat, double val) {
    if ((dist_mat == NULL) || (star1 >= N_stars) || (star2 >= N_stars)) {
        return incorrect_input;
    }
    if (star1 == star2) {
        return ok;
    }
    if (star1 > star2) {
        int tmp = star1;
        star1 = star2;
        star2 = tmp;
    }
    dist_mat[(int)(star1 * (N_stars - 1.5 - 0.5 * star1)) + star2 - 1] = val;
    return ok;
}


star_tracker_error get_dist(int star1, int star2, int N_stars, const double const *dist_mat, double *res) {
    if ((dist_mat == NULL) || (res == NULL) || (star1 >= N_stars) || (star2 >= N_stars)) {
        return incorrect_input;
    }
    if (star1 == star2) {
        *res = 0;
        return ok;
    }
    if (star1 > star2) {
        int tmp = star1;
        star1 = star2;
        star2 = tmp;
    }
    *res = dist_mat[(int)(star1 * (N_stars - 1.5 - 0.5 * star1)) + star2 - 1];
    return ok;
}




// ####################################


void read_image(char *filename, unsigned int *width, unsigned int *height, unsigned char **image) {
    unsigned int error;
    error = lodepng_decode32_file(image, width, height, filename);
    if (error) printf("decoder error %u: %s\n", error, lodepng_error_text(error));
}


void image_to_mat(int w, int h, unsigned char *image, int **mat) {
    int i, j, r, g, b;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            r = image[(i * w + j) * 4];
            g = image[(i * w + j) * 4 + 1];
            b = image[(i * w + j) * 4 + 2];
            mat[i][j] = (int)(0.2126 * r + 0.7152 * g + 0.0722 * b);
        }
    }
}


void mat_to_image(int w, int h, int **mat, unsigned char *image) {
    int i, j;
    //unsigned char *image = malloc(4 * w * h * sizeof(unsigned char));
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            image[(i * w + j) * 4] = mat[i][j];
            image[(i * w + j) * 4 + 1] = mat[i][j];
            image[(i * w + j) * 4 + 2] = mat[i][j];
            image[(i * w + j) * 4 + 3] = 255;
        }
    }
}


void write_image(char *filename, unsigned int width, unsigned int height, unsigned char *image) {
    lodepng_encode32_file(filename, image, width, height);
    //free(image);
}


void gauss_filter(int width, int height, int **pic) {
    int res[height][width], i, j;
    for (i = 1; i < height - 1; i++) {
        for (j = 1; j < width - 1; j++) {
            res[i][j] = 0.084 * pic[i][j]
                      + 0.084 * pic[i + 1][j] + 0.084 * pic[i - 1][j]
                      + 0.084 * pic[i][j + 1] + 0.084 * pic[i][j - 1]
                      + 0.063 * pic[i - 1][j - 1] + 0.063 * pic[i - 1][j + 1]
                      + 0.063 * pic[i + 1][j - 1] + 0.063 * pic[i + 1][j + 1];
        }
    }
    for (i = 1; i < height - 1; i++) {
        for (j = 1; j < width - 1; j++) {
            pic[i][j] = res[i][j];
        }
    }
}


void sobel_filter(int width, int height, int **pic) {
    int res[height][width], i, j, gx, gy;
    for (i = 1; i < height - 1; i++) {
        for (j = 1; j < width - 1; j++) {
            gx = -pic[i - 1][j - 1] + pic[i - 1][j + 1]
                - 2 * pic[i][j - 1] + 2 * pic[i][j + 1]
                - pic[i + 1][j - 1] + pic[i + 1][j + 1];
            gy = pic[i - 1][j - 1] + pic[i - 1][j + 1]
                + 2 * pic[i - 1][j] - 2 * pic[i + 1][j]
                - pic[i + 1][j - 1] - pic[i + 1][j + 1];
            res[i][j] = sqrt(gx * gx + gy * gy);
        }
    }
    for (i = 1; i < height - 1; i++) {
        for (j = 1; j < width - 1; j++) {
            pic[i][j] = res[i][j];
        }
    }
}


star_tracker_error component_search(int width, int height, int **image, int x_start, int y_start, int min_pixel_brightness,
                                    int *comp_size, int *x_mean, int *y_mean, int *sum_brightness) {

    if ((image == NULL) || (comp_size == NULL) || (x_mean == NULL) || (y_mean == NULL) || (sum_brightness == NULL)) {
        return incorrect_input;
    }

    int i = y_start, j, x_left = x_start, x_right = x_start, f = 1;
    *comp_size = 0;
    *sum_brightness = 0;
    *x_mean = 0;
    *y_mean = 0;
    
    while ((i < height - 1) && f) {
        for (j = x_start; (j >= 0) && (image[i][j] > min_pixel_brightness); j--);
        x_left = j + 1;
        for (j = x_start + 1; (j < width) && (image[i][j] > min_pixel_brightness); j++);
        x_right = j - 1;

        f = 0;
        for (j = x_left; j <= x_right; j++) {
            if (image[i][j] > min_pixel_brightness) {
                (*comp_size)++;
                *sum_brightness += image[i][j];
                *x_mean += j * image[i][j];
                *y_mean += i * image[i][j];
                image[i][j] = 0;
            }
            if (image[i + 1][j] > min_pixel_brightness) {
                f = 1;
            }
        }
        i++;
    }

    if (*comp_size) {
        *x_mean /= *sum_brightness;
        *y_mean /= *sum_brightness;
    }    
    return ok;
}


star_tracker_error add_star(list_node **x, list_node **y, list_node **brightness, int x_c, int y_c, int b) {
    if ((x == NULL) || (y == NULL) || (brightness == NULL)) {
        return incorrect_input;
    }
    if ((*brightness == NULL) || (b > (*brightness)->val)) {
        star_tracker_error err;
        if ((err = list_insert(x, x_c)) != ok) {
            return err;
        }
        if ((err = list_insert(y, y_c)) != ok) {
            return err;
        }
        if ((err = list_insert(brightness, b)) != ok) {
            return err;
        }
        return ok;
    }
    return add_star(&((*x)->next), &((*y)->next), &((*brightness)->next), x_c, y_c, b);
}


star_tracker_error detect_stars(int width, int height, int **image, int min_star_size, int max_star_size,
                                int min_pixel_brightness, int min_star_brightness, int max_star_brightness,
                                list_node **x, list_node **y, list_node **brightness) {
    
    if ((image == NULL) || (x == NULL) || (y == NULL) || (brightness == NULL)) {
        return incorrect_input;
    }

    int i, j, size, x_c, y_c, b, n_stars = 0;
    star_tracker_error err;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            if (image[i][j] > min_pixel_brightness) {
                if ((err = component_search(width, height, image, j, i, min_pixel_brightness, &size, &x_c, &y_c, &b)) != ok) {
                    return err;
                }
                if ((size >= min_star_size) && (size <= max_star_size) && (b >= min_star_brightness) && (b <= max_star_brightness)) {
                    n_stars++;
                    if ((err = add_star(x, y, brightness, x_c, y_c, b)) != ok) {
                        return err;
                    }
                }
            }
        }
    }
    if (n_stars < 4) {
        return no_stars_detected;
    }
    return ok;
}


void mark_stars(int w, int h, int **image, list_node *x, list_node *y) {
    int i, j, k;
    for (i = 0; i < CONFIG_SIZE; i++) {

        for (j = -5; j < 5; j++) {
            for (k = -5; k < 5; k++) {
                if ((x->val + j >= 0) && (x->val + j < w) && (y->val + k >= 0) && (y->val + k < h)) {
                    image[y->val + k][x->val + j] = 255;
                }
                
            }
        }

        x = x->next;
        y = y->next;
    }
}


double det_3x3(double a11, double a12, double a13,
               double a21, double a22, double a23,
               double a31, double a32, double a33) {
    return a11 * (a22 * a33 - a32 * a23) - a12 * (a21 * a33 - a31 * a23) + a13 * (a21 * a32 - a31 * a22);
}


star_tracker_error kramer_solve_3x3(double a11, double a12, double a13, double b1,
                                    double a21, double a22, double a23, double b2,
                                    double a31, double a32, double a33, double b3, double *x1, double *x2, double *x3) {

    if ((x1 == NULL) || (x2 == NULL) || (x3 == NULL)) {
        return incorrect_input;
    }

    double delta = det_3x3(a11, a12, a13,
                           a21, a22, a23,
                           a31, a32, a33);
    double delta1 = det_3x3(b1, a12, a13,
                            b2, a22, a23,
                            b3, a32, a33);
    double delta2 = det_3x3(a11, b1, a13,
                            a21, b2, a23,
                            a31, b3, a33);
    double delta3 = det_3x3(a11, a12, b1,
                            a21, a22, b2,
                            a31, a32, b3);
    if (delta == 0) {
        return incorrect_input;
    }
    *x1 = delta1 / delta;
    *x2 = delta2 / delta;
    *x3 = delta3 / delta;
    return ok;
}


star_tracker_error center_coordinates(double *alpha, double *delta, double *r, double *alpha_c, double *delta_c) {
    if ((alpha == NULL) || (delta == NULL) || (r == NULL) || (alpha_c == NULL) || (delta_c == NULL)) {
        return incorrect_input;
    }
    int i;
    double A[4], B[4], C[4], D[4], x[4], y[4], z[4];
    for (i = 0; i < 4; i++) {
        A[i] = cos(delta[i]) * cos(alpha[i]);
        B[i] = cos(delta[i]) * sin(alpha[i]);
        C[i] = sin(delta[i]);
        D[i] = -cos(r[i]);
    }

    kramer_solve_3x3(A[1], B[1], C[1], -D[1],
                     A[2], B[2], C[2], -D[2],
                     A[3], B[3], C[3], -D[3], x + 0, y + 0, z + 0);
    kramer_solve_3x3(A[0], B[0], C[0], -D[0],
                     A[2], B[2], C[2], -D[2],
                     A[3], B[3], C[3], -D[3], x + 1, y + 1, z + 1);
    kramer_solve_3x3(A[0], B[0], C[0], -D[0],
                     A[1], B[1], C[1], -D[1],
                     A[3], B[3], C[3], -D[3], x + 2, y + 2, z + 2);
    kramer_solve_3x3(A[0], B[0], C[0], -D[0],
                     A[1], B[1], C[1], -D[1],
                     A[2], B[2], C[2], -D[2], x + 3, y + 3, z + 3);
    
    double x_mean = 0, y_mean = 0, z_mean = 0;
    for (i = 0; i < 4; i++) {
        x_mean += x[i];
        y_mean += y[i];
        z_mean += z[i];
    }
    double norm = sqrt(x_mean * x_mean + y_mean * y_mean + z_mean * z_mean);
    x_mean /= norm;
    y_mean /= norm;
    z_mean /= norm;

    *delta_c = asin(z_mean);
    *alpha_c = atan2(y_mean, x_mean);
    if (*alpha_c < 0) {
        *alpha_c += 2 * PI;
    }
    return ok;
}


star_tracker_error rotation_angle(double *alpha, double *delta, double *theta, double alpha_c, double delta_c, double *res) {
    if ((alpha == NULL) || (delta == NULL) || (res == NULL)) {
        return incorrect_input;
    }
    double x, y, z, x_b, y_b, x_c, y_c, z_c, n1x, n1y, n1z, n2x, n2y, n2z, theta_true, mx, my, mz;
    int i;

    x_b = cos(alpha_c - 0.5 * PI);
    y_b = sin(alpha_c - 0.5 * PI);

    x_c = cos(delta_c) * cos(alpha_c);
    y_c = cos(delta_c) * sin(alpha_c);
    z_c = sin(delta_c);

    n2x = y_b * z_c;
    n2y = -x_b * z_c;
    n2z = x_b * y_c - y_b * x_c;

    *res = 0;
    double dt;
    for (i = 0; i < 4; i++) {
        x = cos(delta[i]) * cos(alpha[i]);
        y = cos(delta[i]) * sin(alpha[i]);
        z = sin(delta[i]);

        n1x = y * z_c - z * y_c;
        n1y = z * x_c - x * z_c;
        n1z = x * y_c - y * x_c;

        theta_true = acos((n1x * n2x + n1y * n2y + n1z * n2z)
                                / sqrt((n1x * n1x + n1y * n1y + n1z * n1z) * (n2x * n2x + n2y * n2y + n2z * n2z)));
        mx = n1y * n2z - n1z * n2y;
        my = n1z * n2x - n1x * n2z;
        mz = n1x * n2y - n1y * n2x;

        if (mx * x_c + my * y_c + mz * z_c < 0) {
            theta_true = 2 * PI - theta_true;
        }
        dt = theta[i] - theta_true;
        if (dt < 0) {
            dt += 2 * PI;
        }

        *res += dt;
    }

    *res /= 4;
    return ok;
}


star_tracker_error orientation(configuration C, int x_c, int y_c, double deg_per_pixel, int *etalone_configuration, const double const *RA, const double const *dec) {
    double alpha[CONFIG_SIZE], delta[CONFIG_SIZE], r[CONFIG_SIZE], alpha_c, delta_c;
    int i;
    for (i = 0; i < 4; i++) {
        alpha[i] = RA[etalone_configuration[i]] * PI / 180;
        delta[i] = dec[etalone_configuration[i]] * PI / 180;
    }
    for (i = 0; i < CONFIG_SIZE; i++) {
        r[i] = sqrt((C.x[i] - x_c) * (C.x[i] - x_c) + (C.y[i] - y_c) * (C.y[i] - y_c)) * deg_per_pixel * PI / 180;
    }

    center_coordinates(alpha, delta, r, &alpha_c, &delta_c);

    double theta[CONFIG_SIZE], rot;
    for (i = 0; i < CONFIG_SIZE; i++) {
        theta[i] = atan2(y_c - C.y[i], C.x[i] - x_c);
        if (theta[i] < 0) {
            theta[i] += 2 * PI;
        }
    }

    rotation_angle(alpha, delta, theta, alpha_c, delta_c, &rot);

    printf("%lf %lf %lf\n", alpha_c * 180 / PI, delta_c * 180 / PI, rot * 180 / PI);
}



// ####################################


// ####################################


star_tracker_error form_configuration(list_node *x, list_node *y, list_node *b, double deg_per_pixel, configuration *res) {
    if (res == NULL) {
        return incorrect_input;
    }

    star_tracker_error err;
    int i, j;
    for (i = 0; i < CONFIG_SIZE; i++) {
        if ((x == NULL) || (y == NULL) || (b == NULL)) {
            return incorrect_input;
        }

        res->x[i] = x->val;
        res->y[i] = y->val;
        res->b[i] = b->val;

        x = x->next;
        y = y->next;
        b = b->next;
    }

    for (i = 0; i < CONFIG_SIZE; i++) {
        for (j = i + 1; j < CONFIG_SIZE; j++) {
            err = set_dist(i, j, CONFIG_SIZE, res->dist, sqrt((res->x[i] - res->x[j]) * (res->x[i] - res->x[j]) +
                                                        (res->y[i] - res->y[j]) * (res->y[i] - res->y[j])) * deg_per_pixel);
            if (err != ok) {
                return err;
            }
        }
    }

    return ok;
}


star_tracker_error etalone_pairs(double dist, double dist_error, int *first_index, int *pairs_number) {
    if ((first_index == NULL) || (pairs_number == NULL)) {
        return incorrect_input;
    }

    int h1 = hash_table[(int)((dist - 0.5 * dist_error) / DELTA_HASH)];
    while ((h1 < N_PAIRS) && (d[h1] < dist - 0.5 * dist_error)) {
        h1++;
    }
    int h2 = hash_table[(int)((dist + 0.5 * dist_error) / DELTA_HASH)];
    while ((h2 < N_PAIRS) && (d[h1] > dist + 0.5 * dist_error)) {
        h2--;
    }

    *first_index = h1;
    *pairs_number = h2 - h1 + 1;
    return ok;
}


star_tracker_error form_base_table(configuration C, double dist_error) {
    int i, j, first_index, pairs_number, pair;
    star_tracker_error err;
    for (i = 0; i < N_STARS; i++) {
        for (j = 0; j < CONFIG_SIZE - 1; j++) {
            base_table[i][j] = NULL;
        }
    }

    for (i = 0; i < CONFIG_SIZE - 1; i++) {
        if ((err = etalone_pairs(C.dist[i], dist_error, &first_index, &pairs_number)) != ok) {
            return err;
        }
        for (j = first_index; j <= first_index + pairs_number; j++) {
            if ((err = list_insert(base_table[k[j]] + i, l[j])) != ok) {
                return err;
            }
            if ((err = list_insert(base_table[l[j]] + i, k[j])) != ok) {
                return err;
            }
            // if ((k[j] == 739) || (l[j] == 739)) {
            //     printf("%d %d %lf\n", k[j], l[j], d[j]);
            // }
        }
    }
    return ok;
}


int all_not_null(list_node **ptr_arr) {
    int i;
    for (i = 0; i < CONFIG_SIZE; i++) {
        if (ptr_arr[i] == NULL) {
            return 0;
        }
    }
    return 1;
}


star_tracker_error increment_ptr_arr(list_node **base_table_line, list_node **ptr_arr) {
    if ((base_table_line == NULL) || (ptr_arr == NULL)) {
        return incorrect_input;
    }
    int i = 1;
    ptr_arr[i] = ptr_arr[i]->next;
    while ((i < CONFIG_SIZE - 1) && (ptr_arr[i] == NULL)) {
        ptr_arr[i] = base_table_line[i - 1];
        ptr_arr[i + 1] = ptr_arr[i + 1]->next;
        i++;
    }
    if (ptr_arr[i] == NULL) {
        return incorrect_input;
    }
    return ok;
}


star_tracker_error add_configuration_line(int ***star_num_table, double ***star_dist_table, int *table_size, int *curr_size,
                                          list_node **ptr_arr) {
    int i, j;
    star_tracker_error err;
    int **tmp1, m;
    double **tmp2, dist;

    if (*table_size >= *curr_size) {
        *curr_size *= 2;
        tmp1 = realloc(*star_num_table, *curr_size * sizeof(int *));
        if (tmp1 == NULL) {
            return memory_error;
        }
        *star_num_table = tmp1;
        for (j = *table_size; j < *curr_size; j++) {
            (*star_num_table)[j] = NULL;
        }

        tmp2 = realloc(*star_dist_table, *curr_size * sizeof(double *));
        if (tmp2 == NULL) {
            return memory_error;
        }
        *star_dist_table = tmp2;
        for (j = *table_size; j < *curr_size; j++) {
            (*star_dist_table)[j] = NULL;
        }
    }

    (*star_num_table)[*table_size] = malloc(CONFIG_SIZE * sizeof(int));
    if ((*star_num_table)[*table_size] == NULL) {
        return memory_error;
    }
    for (j = 0; j < CONFIG_SIZE; j++) {
        (*star_num_table)[*table_size][j] = ptr_arr[j]->val;
    }
    (*star_dist_table)[*table_size] = malloc((CONFIG_SIZE * (CONFIG_SIZE - 1) / 2) * sizeof(double));
    if ((*star_dist_table)[*table_size] == NULL) {
        return memory_error;
    }

    for (j = 0; j < CONFIG_SIZE; j++) {
        for (m = j + 1; m < CONFIG_SIZE; m++) {
            //printf("%d %d\n", j, m);
            err = get_dist(ptr_arr[j]->val, ptr_arr[m]->val, N_STARS, dist_mat, &dist);
            if (err != ok) {
                return err;
            }
            //printf("e\n");
            err = set_dist(j, m, CONFIG_SIZE, (*star_dist_table)[*table_size], dist);
            if (err != ok) {
                return err;
            }
        }
    }
    //printf("f\n");

    (*table_size)++;
    //printf("g\n");
    return ok;
}


star_tracker_error form_configuration_table(int ***star_num_table, double ***star_dist_table, int *table_size) {
    int i, j, curr_size = 1;
    *table_size = 0;
    if ((star_num_table == NULL) || (star_dist_table == NULL) || (table_size == NULL)) {
        return incorrect_input;
    }
    *star_num_table = calloc(curr_size, sizeof(int *));
    *star_dist_table = calloc(curr_size, sizeof(double *));
    if ((star_dist_table == NULL) || (star_num_table == NULL)) {
        return memory_error;
    }

    star_tracker_error err;
    list_node *ptr_arr[CONFIG_SIZE], ln0, *tmp;
    ln0.next = NULL;
    for (i = 0; i < N_STARS; i++) {
        if (all_not_null(base_table[i])) {
            ptr_arr[0] = &ln0;
            ptr_arr[0]->val = i;
            for (j = 0; j < CONFIG_SIZE - 1; j++) {
                ptr_arr[j + 1] = base_table[i][j];
            }
            
            add_configuration_line(star_num_table, star_dist_table, table_size, &curr_size, ptr_arr);
        

            while (increment_ptr_arr(base_table[i], ptr_arr) == ok) {
                add_configuration_line(star_num_table, star_dist_table, table_size, &curr_size, ptr_arr);
            }
        }
    }

    return ok;
}


star_tracker_error select_configuration(configuration C, int table_size, int **star_num_table, double **star_dist_table,
                                        double max_sq_err, list_node **res, int *res_num) {

    if ((star_num_table == NULL) || (star_dist_table == NULL) || (res == NULL) || (res_num == NULL)) {
        return incorrect_input;
    }
    int i, j;
    double sq_err;
    *res_num = 0;
    star_tracker_error err;
    for (i = 0; i < table_size; i++) {
        sq_err = 0;
        for (j = 0; j < CONFIG_SIZE * (CONFIG_SIZE - 1) / 2; j++) {
            sq_err += (star_dist_table[i][j] - C.dist[j]) * (star_dist_table[i][j] - C.dist[j]);
        }
        if ((sq_err <= max_sq_err) /*&& (mag[star_num_table[i][0]] <= mag[star_num_table[i][1]])
                                   && (mag[star_num_table[i][0]] <= mag[star_num_table[i][2]])
                                   && (mag[star_num_table[i][0]] <= mag[star_num_table[i][3]])*/) {
            printf("%d %d %d %d\n", star_num_table[i][0], star_num_table[i][1], star_num_table[i][2], star_num_table[i][3]);
            printf("%lf\n", sq_err);
            if ((err = list_insert(res, i)) != ok) {
                return err;
            }
            (*res_num)++;
        }
    }
    return ok;
} 

