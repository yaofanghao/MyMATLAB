/*
* Copyright (c) 2016, Gabor Adam Fodor <fogggab@yahoo.com>
* All rights reserved.
*
* License:
*
* This program is provided for scientific and educational purposed only.
* Feel free to use and/or modify it for such purposes, but you are kindly
* asked not to redistribute this or derivative works in source or executable
* form. A license must be obtained from the author of the code for any other use.
*
*/
#include "FFST.h"
#include <iostream>
#include <fstream>

/**
* Fast finite shearlet tranform implementation following the Matlab implementation
* of the image processing group of TU Kaiserslautern. The matlab code and detailed
* publications on their method can be found at:
* http://www.mathematik.uni-kl.de/imagepro/software/ffst/
* A copy of the same matlab files is also available via the Platypus project at
* http://www.project-platypus.net/cradledownload.html
*
* NOTE: This implementation ONLY works for images of size 512x512, skipping over
* generating the filter-banks on the run to speed up computation speed. Instead,
* all frequency domain filters are stored in a big lookup table that is initalized
* the first time the transform is called.
*
* For a more detailed explanation on the functions in this file, check out the Matlab
* implementation.
**/

namespace FFST{

	std::vector<cv::Mat> Psi;

	cv::Mat inverseShearletTransformSpect(std::vector<cv::Mat> &ST){

		//Only supports 512x512 images
		cv::Mat img(512, 512, CV_32F, cv::Scalar(0));
		
		//Check if Psi is initialized
		if (Psi.size() == 0){
			//Read coefficients
			set_coeffs_512x512();
		}

		cv::Mat dftout, scc, src;
		scc = cv::Mat(512, 512, CV_32F, cv::Scalar(0));
		src = cv::Mat(512, 512, CV_32F, cv::Scalar(0));

		for (int k = 0; k < ST.size(); k++){
			cv::Mat cc, rc;
			cv::dft(ST[k], dftout, cv::DFT_COMPLEX_OUTPUT);

			cc = cv::Mat(512, 512, CV_32F, cv::Scalar(0));
			rc = cv::Mat(512, 512, CV_32F, cv::Scalar(0));

			//Separate complex output
			for (int i = 0; i < cc.rows; i++){
				for (int j = 0; j < cc.cols; j++){
					rc.at<float>(i, j) = dftout.at<float>(i, j * 2);
					cc.at<float>(i, j) = dftout.at<float>(i, j * 2 + 1);
				}
			}

			rc = fftshift(rc);
			cc = fftshift(cc);

			for (int i = 0; i < cc.rows; i++){
				for (int j = 0; j < cc.cols; j++){
					rc.at<float>(i, j) *= Psi[k].at<float>(i, j);
					cc.at<float>(i, j) *= Psi[k].at<float>(i, j);
				}
			}

			src += rc;
			scc += cc;
		}
		
		src = ifftshift(src);
		scc = ifftshift(scc);

		//Put it all together for forward transform
		cv::Mat xx(src.rows, src.cols, CV_32FC2, cv::Scalar(0));
		for (int i = 0; i < xx.rows; i++){
			for (int j = 0; j < xx.cols; j++){
				xx.at<cv::Point2f>(i, j).x = src.at<float>(i, j);	//Real part
				xx.at<cv::Point2f>(i, j).y = scc.at<float>(i, j);	//Imaginary part
			}
		}

		cv::dft(xx, img, cv::DFT_SCALE | cv::DFT_INVERSE | cv::DFT_REAL_OUTPUT);

		return img;
	}

	std::vector<cv::Mat> shearletTransformSpect(cv::Mat &img){
		std::vector<cv::Mat> ST(61);

		if (img.rows != 512 || img.cols != 512){
			std::cerr << "Image resolution not supported!" << std::endl;
		}

		//Check if Psi is initialized
		if (Psi.size() == 0){
			//Read coefficients
			set_coeffs_512x512();
		}

		//fftshift(fft2(A))
		cv::Mat dftout, cc, rc;
		cv::dft(img, dftout, cv::DFT_COMPLEX_OUTPUT);

		//Separate complex and real channel
		cc = cv::Mat(img.rows, img.cols, CV_32F, cv::Scalar(0));
		rc = cv::Mat(img.rows, img.cols, CV_32F, cv::Scalar(0));

		for (int i = 0; i < img.rows; i++){
			for (int j = 0; j < img.cols; j++){
				rc.at<float>(i, j) = dftout.at<float>(i, j * 2);
				cc.at<float>(i, j) = dftout.at<float>(i, j * 2 + 1);
			}
		}
		rc = fftshift(rc);
		cc = fftshift(cc);

		for (int k = 0; k < 61; k++){
			cv::Mat lrc, lcc;
			lcc = cv::Mat(cc.rows, cc.cols, CV_32F);
			lrc = cv::Mat(rc.rows, rc.cols, CV_32F);

			for (int i = 0; i < rc.rows; i++){
				for (int j = 0; j < rc.cols; j++){
					// Psi(:,:,k) .* fft2(img)
					lrc.at<float>(i, j) = Psi[k].at<float>(i, j) * rc.at<float>(i, j);
					lcc.at<float>(i, j) = Psi[k].at<float>(i, j) * cc.at<float>(i, j);
				}
			}

			lrc = fftshift(lrc);
			lcc = fftshift(lcc);

			//Put it all together for forward transform
			cv::Mat xx(lrc.rows, lrc.cols, CV_32FC2, cv::Scalar(0));
			for (int i = 0; i < xx.rows; i++){
				for (int j = 0; j < xx.cols; j++){
					xx.at<cv::Point2f>(i, j).x = lrc.at<float>(i, j);	//Real part
					xx.at<cv::Point2f>(i, j).y = lcc.at<float>(i, j);	//Imaginary part
				}
			}

			cv::dft(xx, ST[k], cv::DFT_SCALE | cv::DFT_INVERSE | cv::DFT_REAL_OUTPUT);
		}

		return ST;
	}

	void set_coeffs_512x512(){

		std::ifstream myfile;
		myfile.open("FFST_512x512_table.txt");

		if (!myfile){
			std::cerr << "Cannot find coefficient file FFST_512x512_table.txt! " << std::endl;
			return;
		}

		Psi = std::vector<cv::Mat>(61);
		for (int i = 0; i < Psi.size(); i++){
			Psi[i] = cv::Mat(512, 512, CV_32F, cv::Scalar(0));
		}
		int n;
		myfile >> n; //Number of lines to read
		for (int i = 0; i<n; i++){
			int a, b, c;
			double v;

			myfile >> a >> b >> c >> v;
			Psi[c - 1].at<float>(a - 1, b - 1) = v;
		}

		myfile.close();
	}
	
	cv::Mat fftshift(cv::Mat &x){
		//Inverse FFT shift
		// FFTSHIFT(X) swaps the first and third
		// quadrants and the second and fourth quadrants.For N - D
		// arrays, FFTSHIFT(X) swaps "half-spaces" of X along each
		// dimension.

		return circshift(x, std::floor(x.rows / 2), std::floor(x.cols / 2));
	}

	cv::Mat ifftshift(cv::Mat &x){
		//Inverse FFT shift
		// IFFTSHIFT(X) swaps the first and third
		// quadrants and the second and fourth quadrants.For N - D
		// arrays, IFFTSHIFT(X) swaps "half-spaces" of X along each
		// dimension.

		return circshift(x, std::ceil(x.rows * 1.0 / 2), std::ceil(x.cols * 1.0 / 2));
	}

	cv::Mat circshift(cv::Mat &M, int x, int y){
		// Circularly shifts the values in matrix M
		// by x along the first dimension and by y along the second one

		if ((x == 0) && (y == 0))
			return M;

		int X = M.rows;
		int Y = M.cols;

		cv::Mat res(X, Y, CV_32F, cv::Scalar(0));
		for (int i = 0; i < M.rows; i++){
			for (int j = 0; j < M.cols; j++){
				res.at<float>((i + x + X) % X, (j + y + Y) % Y) = M.at<float>(i, j);
			}
		}
		return res;
	}
}