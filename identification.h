#ifndef IDENTIFICATION_H
#define IDENTIFICATION_H


#include "catalogs.h"

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "lodepng.h"


typedef enum star_tracker_error {
    ok,
    incorrect_input,
    no_stars_detected,
    memory_error,
    zero_configurations,
    multiple_configurations
} star_tracker_error;


typedef struct list_node{
    int val;
    struct list_node* next;
} list_node;


typedef struct configuration {
    int x[CONFIG_SIZE], y[CONFIG_SIZE], b[CONFIG_SIZE];
    double dist[CONFIG_SIZE * (CONFIG_SIZE - 1) / 2];
} configuration;


star_tracker_error list_insert(list_node **first, int val);

star_tracker_error list_clear(list_node **first);

extern list_node* base_table[N_STARS][CONFIG_SIZE - 1];

star_tracker_error set_dist(int star1, int star2, int N_stars, double *dist_mat, double val);

star_tracker_error get_dist(int star1, int star2, int N_stars, const double const *dist_mat, double *res);






void read_image(char *filename, unsigned int *width, unsigned int *height, unsigned char **image);

void image_to_mat(int w, int h, unsigned char *image, int **mat);

void mat_to_image(int w, int h, int **mat, unsigned char *image);

void write_image(char *filename, unsigned int width, unsigned int height, unsigned char *image);

void gauss_filter(int width, int height, int **pic);

void sobel_filter(int width, int height, int **pic);

star_tracker_error component_search(int width, int height, int **image, int x_start, int y_start, int min_pixel_brightness,
    int *comp_size, int *x_mean, int *y_mean, int *sum_brightness);

star_tracker_error add_star(list_node **x, list_node **y, list_node **brightness, int x_c, int y_c, int b);

star_tracker_error detect_stars(int width, int height, int **image, int min_star_size, int max_star_size,
    int min_pixel_brightness, int min_star_brightness, int max_star_brightness,
    list_node **x, list_node **y, list_node **brightness);

void mark_stars(int w, int h, int **image, list_node *x, list_node *y);

double det_3x3(double a11, double a12, double a13,
               double a21, double a22, double a23,
               double a31, double a32, double a33);

star_tracker_error kramer_solve_3x3(double a11, double a12, double a13, double b1,
                                    double a21, double a22, double a23, double b2,
                                    double a31, double a32, double a33, double b3, double *x1, double *x2, double *x3);

star_tracker_error center_coordinates(double *alpha, double *delta, double *r, double *alpha_c, double *delta_c);

star_tracker_error rotation_angle(double *alpha, double *delta, double *theta, double alpha_c, double delta_c, double *res);

star_tracker_error orientation(configuration C, int x_c, int y_c, double deg_per_pixel, int *etalone_configuration,
                               const double const *RA, const double const *dec);





// ##########################

// ##########################


star_tracker_error form_configuration(list_node *x, list_node *y, list_node *b, double deg_per_pixel, configuration *res);

star_tracker_error etalone_pairs(double dist, double dist_error, int *first_index, int *pairs_number);

star_tracker_error form_base_table(configuration C, double dist_error);

int all_not_null(list_node **ptr_arr);

star_tracker_error increment_ptr_arr(list_node **base_table_line, list_node **ptr_arr);

star_tracker_error add_configuration_line(int ***star_num_table, double ***star_dist_table, int *table_size, int *curr_size,
                                          list_node **ptr_arr);

star_tracker_error form_configuration_table(int ***star_num_table, double ***star_dist_table, int *table_size);

star_tracker_error select_configuration(configuration C, int table_size, int **star_num_table, double **star_dist_table,
                                        double max_sq_err, list_node **res, int *res_num);

#endif