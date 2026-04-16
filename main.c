#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "identification.h"


void clear_memory(int height, int **image, list_node **x, list_node **y, list_node **b,
                  int table_size, int **star_num_table, double **star_dist_table) {
    int i, j;

    if (image != NULL) {
        for (i = 0; i < height; i++) {
            if (image[i] != NULL) {
                free(image[i]);
            }
        }
        free(image);
    }

    if (x != NULL) {
        list_clear(x);
    }
    if (y != NULL) {
        list_clear(y);
    }
    if (b != NULL) {
        list_clear(b);
    }

    for (i = 0; i < N_STARS; i++) {
        for (j = 0; j < CONFIG_SIZE - 1; j++) {
            list_clear(base_table[i] + j);
        }
    }

    if (star_num_table != NULL) {
        for (i = 0; i < table_size; i++) {
            if (star_num_table[i] != NULL) {
                free(star_num_table[i]);
            }
        }
        free(star_num_table);
    }
    if (star_dist_table != NULL) {
        for (i = 0; i < table_size; i++) {
            if (star_dist_table[i] != NULL) {
                free(star_dist_table[i]);
            }
        }
        free(star_dist_table);
    }
}


int main(int argc, char **argv) {
    char filename[50] = "test_9.png";
    if (argc >= 2) {
        strcpy(filename, argv[1]);
    }

    int i, j;
    unsigned int width, height;
    unsigned char *image_read;
    read_image(filename, &width, &height, &image_read);
    
    int **image = malloc(height * sizeof(int *));
    for (i = 0; i < height; i++) {
        image[i] = malloc(width * sizeof(int));
    }
    image_to_mat(width, height, image_read, image);

    //free(image_read);

    // нижележащий код не зависит от происхождения image

    star_tracker_error err;
    list_node *x = NULL, *y = NULL, *brightness = NULL;
    if ((err = detect_stars(width, height, image, MIN_STAR_SIZE, MAX_STAR_SIZE, MIN_PIXEL_BRIGHTNESS, MIN_STAR_BRIGHTNESS, MAX_STAR_BRIGHTNESS, &x, &y, &brightness)) != ok) {
        //clear_memory(height, image, &x, &y, &brightness, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, NULL, NULL);
        return err;
    }

    mark_stars(width, height, image, x, y);

    mat_to_image(width, height, image, image_read);
    write_image("out.png", width, height, image_read);
    free(image_read);
    
    configuration C;
    if ((err = form_configuration(x, y, brightness, DEG_PER_PIXEL, &C)) != ok) {  // мой телефон - 0.02, камера Александра - 0.01
        //clear_memory(height, image, &x, &y, &brightness, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, NULL, NULL);
        return err;
    }

    for (i = 0; i < CONFIG_SIZE * (CONFIG_SIZE - 1) / 2; i++) {
        printf("%lf\n", C.dist[i]);
    }

    if ((err = form_base_table(C, DIST_ERROR)) != ok) {
        //clear_memory(height, image, &x, &y, &brightness, dist_mat, hash_table, k, l, d, RA, dec, N_stars, base_table, 0, NULL, NULL);
        return err;
    }


    int **star_num_table, table_size;
    double **star_dist_table;
    if ((err = form_configuration_table(&star_num_table, &star_dist_table, &table_size)) != ok) {
        //clear_memory(height, image, &x, &y, &brightness, dist_mat, hash_table, k, l, d, RA, dec, N_stars, candidate_base_table, table_size, star_num_table, star_dist_table);
        return err;
    }

    printf("%d\n", table_size);

    list_node *conf_ind = NULL;
    int conf_num;
    if ((err = select_configuration(C, table_size, star_num_table, star_dist_table, MAX_SQ_ERR, &conf_ind, &conf_num)) != ok) {
        //clear_memory(height, image, &x, &y, &brightness, dist_mat, hash_table, k, l, d, RA, dec, table_size, star_num_table, star_dist_table);
        return err;
    }

    printf("%d\n", conf_num);
    

    // if (conf_num == 0) {

    // }
    // else if (conf_num > 1) {
    //     printf("%d %d\n", C.x[2], C.y[2]);
    //     update_configuration(x, y, brightness, 0.01, &C);
    //     printf("%d %d\n", C.x[2], C.y[2]);
    //     for (i = 0; i < N_stars; i++) {
    //         if (candidate_base_table[i] != NULL) {
    //             for (j = 0; j < 3; j++) {
    //                 list_clear(candidate_base_table[i] + j);
    //             }
    //             free(candidate_base_table[i]);
    //         }
    //     }
    //     form_candidate_base_table(C, dist_error, max_angle / hash_table_size, hash_table, N_pairs, k, l, d, N_stars, candidate_base_table);

    //     int **star_num_table_1, table_size_1;
    //     double **star_dist_table_1;
    //     form_candidate_configuration_table(N_stars, candidate_base_table, dist_mat, &star_num_table_1, &star_dist_table_1, &table_size_1);

    //     list_node *conf_ind_1 = NULL;
    //     int conf_num_1;
    //     select_configuration(C, table_size_1, star_num_table_1, star_dist_table_1, mag, 0.5, &conf_ind_1, &conf_num_1);

    //     list_node *tmp1 = conf_ind, *tmp2;
    //     while (tmp1 != NULL) {
    //         tmp2 = conf_ind_1;
    //         while (tmp2 != NULL) {
    //             if ((star_num_table[tmp1->val][0] == star_num_table_1[tmp2->val][0])) {
    //                 printf("%d %d %d %d\n", star_num_table[tmp1->val][0], star_num_table[tmp1->val][1], star_num_table[tmp1->val][2], star_num_table[tmp1->val][3]);
    //             }
    //             tmp2 = tmp2->next;
    //         }
    //         tmp1 = tmp1->next;
    //     }
    // }

    //printf("%d\n", configuration_index);

    // double tmp;
    // get_dist(1022, 1057, N_stars, dist_mat, &tmp);
    // printf("%lf\n", tmp);
    

    orientation(C, width / 2, height / 2, DEG_PER_PIXEL, star_num_table[conf_ind->val], RA, dec);//here

    //clear_memory(height, image, &x, &y, &brightness, dist_mat, hash_table, k, l, d, RA, dec, N_stars, candidate_base_table, table_size, star_num_table, star_dist_table);
    
    // list_clear(&conf_ind);

    return 0;
}

//gcc main.c lodepng.c lodepng.h catalogs.h catalogs.c identification.c identification.h -lm -o main