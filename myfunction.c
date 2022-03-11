#include <stdbool.h>
#include <malloc.h>
#include "readBMP.h"
#include "writeBMP.h"


typedef struct {
    unsigned char red;
    unsigned char green;
    unsigned char blue;

    // add in order to make struct size 4 bytes
    unsigned char offset;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;

    // add in order to make struct size 16 bytes
    int sum;
} pixel_sum;


/* Compute min and max of two integers, respectively */
int min(int a, int b) { return (a < b ? a : b); }
int max(int a, int b) { return (a > b ? a : b); }

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

    // divide by kernel's weight
    sum.red = sum.red / kernelScale;
    sum.green = sum.green / kernelScale;
    sum.blue = sum.blue / kernelScale;

    // truncate each pixel's color values to match the range [0,255]
    current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
    current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
    current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {
    sum->red += ((int) p.red) * weight;
    sum->green += ((int) p.green) * weight;
    sum->blue += ((int) p.blue) * weight;
    // sum->num++;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
    int iii = max(i-1, 0);
    int jjj = max(j-1, 0);
    int iiLimit = min(i+1, dim-1);
    int jjLimit = min(j+1, dim-1);

    int ii = iii, jj;
    int currRow, currCol;
    pixel_sum sum = {0,0,0,0};
    pixel current_pixel;
    pixel loop_pixel;

    for(; ii <= iiLimit; ii++) {
        for(jj = jjj; jj <= jjLimit; jj++) {

            int kRow, kCol;

            // compute row index in kernel
            if (ii > i) {
                kRow = 2;
            } else if (ii < i) {
                kRow = 0;
            } else {
                kRow = 1;
            }

            // compute column index in kernel
            if (jj > j) {
                kCol = 2;
            } else if (jj < j) {
                kCol = 0;
            } else {
                kCol = 1;
            }

            // apply kernel on pixel at [ii,jj]
            sum_pixels_by_weight(&sum, src[ii * dim + jj], kernel[kRow][kCol]);
        }
    }

    if (filter) {
        int min_row, min_col, max_row, max_col;
        int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
        int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
        // find min and max coordinates
        for(ii = iii; ii <= iiLimit; ii++) {
            for(jj = jjj; jj <= jjLimit; jj++) {
                // check if smaller than min or higher than max and update
                loop_pixel = src[ii * dim + jj];
                int pixelSum = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
                if (pixelSum < min_intensity) {
                    min_intensity = pixelSum;
                    min_row = ii;
                    min_col = jj;
                }
                if (pixelSum > max_intensity) {
                    max_intensity = pixelSum;
                    max_row = ii;
                    max_col = jj;
                }
            }
        }
        // filter out min and max
        sum_pixels_by_weight(&sum, src[min_row * dim + min_col], -1);
        sum_pixels_by_weight(&sum, src[max_row * dim + max_col], -1);
    }

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
    int kernelSizeDivided = kernelSize >> 1;
    int limit =  dim - kernelSizeDivided;

    int i = kernelSizeDivided, j;
    for (; i < limit; i++) {
        for (j =  kernelSizeDivided ; j < limit; j++) {
            dst[i * dim + j] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
        }
    }
}

void charsToPixels(Image *charsImg, pixel* pixels) {

    int row = 0, col;
    for (; row < m ; row++) {
        int rowNum = row*n;
        int tripleRowNum = (rowNum << 1) + rowNum;


        for (col = 0 ; col < n ; col++) {
            int tripleColumn = (col << 1) + col;
            int pixelLocation = rowNum + col;

            pixels[pixelLocation].red = image->data[tripleRowNum + tripleColumn];
            pixels[pixelLocation].green = image->data[tripleRowNum + tripleColumn + 1];
            pixels[pixelLocation].blue = image->data[tripleRowNum + tripleColumn + 2];
        }
    }
}

void pixelsToChars(pixel* pixels, Image *charsImg) {

    int row = 0, col;
    for (; row < m ; row++) {
        int rowNum = row * n;
        int tripleRowNum = ((row << 1) + row) * n;
        for (col = 0 ; col < n ; col+=3) {
            int tripleCol = (col << 1) + col;
            int tripleSecondCol = (col << 1) + col + 4;
            int tripleThirdCol = (col << 1) + col + 5;
            int tripleSum = tripleRowNum + tripleCol;
            int rowNumColSum = rowNum + col;

            image->data[tripleSum] = pixels[rowNumColSum].red;
            image->data[tripleSum + 1] = pixels[rowNumColSum].green;
            image->data[tripleSum + 2] = pixels[rowNumColSum].blue;

            image->data[tripleSum + 3] = pixels[rowNumColSum + 1].red;
            image->data[tripleSum + 4] = pixels[rowNumColSum + 1].green;
            image->data[tripleSum + 5] = pixels[rowNumColSum + 1].blue;

            image->data[tripleSum + 6] = pixels[rowNumColSum + 2].red;
            image->data[tripleSum + 7] = pixels[rowNumColSum + 2].green;
            image->data[tripleSum + 8] = pixels[rowNumColSum + 2].blue;
        }
        for (; col < n ; col++) {
            int tripleCol = (col << 1) + col;
            int tripleSum = tripleRowNum + tripleCol;

            image->data[tripleSum] = pixels[row*n + col].red;
            image->data[tripleSum + 1] = pixels[row*n + col].green;
            image->data[tripleSum + 2] = pixels[row*n + col].blue;
        }
    }
}

void copyPixels(pixel* src, pixel* dst) {

    int row, col;
    for (row = 0 ; row < m ; row++) {
        int rowNum = row * n;
        for (col = 0 ; col < n ; col++) {
            int pixelLocation = rowNum + col;

            dst[pixelLocation].red = src[pixelLocation].red;
            dst[pixelLocation].green = src[pixelLocation].green;
            dst[pixelLocation].blue = src[pixelLocation].blue;
        }
    }
}

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
    int mallocSize = m*n*sizeof(pixel);

    pixel* pixelsImg = malloc(mallocSize);
    pixel* backupOrg = malloc(mallocSize);

    charsToPixels(image, pixelsImg);
    copyPixels(pixelsImg, backupOrg);

    smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

    pixelsToChars(pixelsImg, image);

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

    /*
    * [1, 1, 1]
    * [1, 1, 1]
    * [1, 1, 1]
    */
    int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    /*
    * [-1, -1, -1]
    * [-1, 9, -1]
    * [-1, -1, -1]
    */
    int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};

    if (flag == '1') {
        // blur image
        doConvolution(image, 3, blurKernel, 9, false);

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);

        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName);
    } else {
        // apply extermum filtered kernel to blur image
        doConvolution(image, 3, blurKernel, 7, true);

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName);
    }
}
