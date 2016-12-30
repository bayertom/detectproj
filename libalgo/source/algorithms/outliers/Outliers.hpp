// Description: Outliers detection using the least squares and M-estimators

// Copyright (c) 2010 - 2014
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.


#ifndef Outliers_HPP
#define Outliers_HPP

template <typename Point1, typename Point2, typename TKey>
void Outliers::findOutliersLTS(const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source, Container <Point1 *> &global_points_dest,
	Container <Point2 *> &local_points_dest, TKey & min_key, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, const typename Point1::Type perc_ratio)
{
	//Removing outliers using the iteratively LTS
	//Outliers are dectected using weighted transformation and Danish method

	const unsigned int n_global_points_source = global_points_source.size();
	unsigned int n_iterations = 10 * n_global_points_source;
	bool best_key_longer = false;

	//Throw exception: not enough points
	if (n_global_points_source < 3)
	{
		throw BadDataException("BadDataException: not enough global points to be optimized. ", "Can not find optimal transformation key.");
	}

	//Initialize min error
	typename Point1::Type min_error = MAX_FLOAT;

	//Initialize random number generator
	srand((unsigned)time(0));

	//Perform n iterations to find the best key
	unsigned int i = 0, best_key_size = 2, pairs_actual_size = 0;
	TAccuracyCharacteristics <typename Point1::Type> accuracy;

	do
	{
		//Get 2 randomly generated different indices
		int index2 = 0;
		int index1 = rand() % n_global_points_source;

		do
		{
			index2 = rand() % n_global_points_source;
		} while (index1 == index2);

		//Create point pairs
		typename TDevIndexPairs <typename Point1::Type>::Type pairs, pairs_actual;
		pairs_actual.push_back(std::make_pair(0.0, index1));
		pairs_actual.push_back(std::make_pair(1.0, index2));

		//Rearrange points: use only fist k-best points
		rearrangePoints(global_points_source, local_points_source, global_points_dest, local_points_dest, pairs_actual);

		//Create keys and errors
		TKey key, key_actual;
		typename Point1::Type error = MAX_FLOAT, error_actual = 0.1 * error;

		//Performs iterations for this key until it convergates
		for (unsigned int j = 0; j < 2, error_actual < error; j++)
		{
			//Assign old values
			error = error_actual;
			key = key_actual;
			pairs = pairs_actual;

			//Compute transformation for the key of 2 points
			Container <Point1 *> transformed_points;
			getTransformKey(global_points_dest, local_points_dest, key_actual);
			transform(global_points_source, local_points_source, &transformed_points, key_actual);

			//Compute deviation
			TAccuracyCharacteristics <typename Point1::Type> accuracy;
			accuracy = getAccuracyCharacteristics(global_points_source, local_points_source, &transformed_points, key_actual);

			//Find new k-best pairs and sort
			createKBestPairsOfPoints(accuracy, pairs_actual, perc_ratio);

			//Rearrange points: use only fist k-best points
			rearrangePoints(global_points_source, local_points_source, global_points_dest, local_points_dest, pairs_actual);

			//Compute transformation only k_best local points using k-best key
			Container <Point1 *> transformed_points_dest;
			transform(global_points_dest, local_points_dest, transformed_points_dest, key_actual);

			//Get accuracy characteristic for k-best global points and k-best transformed points using k-best key
			TAccuracyCharacteristics <typename Point1::Type> accuracy_actual2 = getAccuracyCharacteristics(global_points_dest, local_points_dest, &transformed_points_dest, key_actual);
			error_actual = accuracy_actual2.std_dev;
		}

		//We found a better key: remember deviation (error), key and pairs
		if (error_actual < min_error)
		{
			min_error = error_actual;
			min_key = key;
			min_pairs = pairs;
		}

	} while (++i < n_iterations && perc_ratio != 1.0);

	//Sort pairs points indices
	std::sort(min_pairs.begin(), min_pairs.end(), sortPointPairsByIndices <typename Point1::Type>());

	//Rearrange points: use only first k-best points
	rearrangePoints(global_points_source, local_points_source, global_points_dest, local_points_dest, min_pairs);
}


template <typename T>
void Outliers::findOutliersME(const Matrix <T> &P, const Matrix <T> &Q, const T k, const T tol, const TMEstimatorsScheme me_scheme, TMEstimatorsWeightFunction me_weight_function, const unsigned int max_iter, Matrix <T> &W, Matrix <unsigned int> &I, Matrix <T> &Eps,  T &f_init, T &f, unsigned int &iter)
{
	//Removing outliers using the M-estimates: P = correct, Q = contains outliers
	//Solution based on IRLS
	//Scheme and weight functions may be selected
	//Sn based on the MAD (fixed scaling)
	const unsigned int m1 = P.rows(), n1 = P.cols(), m2 = P.rows(), n2 = P.cols(), m3 = W.rows(), n3 = W.cols(), m4 = I.rows(), m5 = Eps.rows();
	const T wmin = 0.001, wmax = 0.9;

	//Test number of rows A and IX
	if (m1 != m2)
	{
		throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different rows count for P, Q, can not find outliers:  ", P, Q);
	}

	//Test number of columns A and IX
	if (n1 != n2)
	{
		throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different columns count for P, Q, can not find outliers:  ", P, Q);
	}

	//W matrix must be square
	if (m3 != n3)
	{
		throw MathMatrixNotSquareException <Matrix <T> >("MatrixMatrixNotSquareException: ", " W matrix must be square:  ", W);
	}

	//Are there enough points?
	if ((me_scheme == ScaleScheme) && (m1 < 2) || (me_scheme == !ScaleScheme) && (m1 < 4))
	{
		throw BadDataException("BadDataException: not enough points, ", " can not find outliers using M-estimators.");
	}

	//W and P
	if (m3 != 2 * m1)
	{
		throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " incorrect size of P(n, n), and W(2*n,2*n), :  ", P, W);

	}

	//W and I
	if (m3 != m4)
	{
		throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different rows count for W, I, can not find outliers:  ", W, I);
	}

	//I and Eps
	if (m4 != m5)
	{
		throw MathMatrixDifferentSizeException <Matrix <T> >("MathMatrixDifferentSizeException: ", " different rows count for I, Eps, can not find outliers:  ", I, Eps);
	}

	//Create A matrix
	unsigned int n_cols = 4;
	if (me_scheme == ShiftsScheme) n_cols = 2;
	else if (me_scheme == ScaleScheme) n_cols = 1;
	else if (me_scheme == ScaleShiftsScheme) n_cols = 3;
	Matrix <T> A(2 * m1, n_cols);

	//Set elements of A(2m, 1) matrix, applicable to each scheme
	Matrix <T> E = ones(m1, 1, 1.0);
	if (me_scheme != ShiftsScheme)
	{
		A.replace(Q(0, m1 - 1, 0, 0), 0, 0);
		A.replace(Q(0, m1 - 1, 1, 1), m1, 0);
	}

	else
	{
		A.replace(E, 0, 0);
		A.replace(E, m1, 1);
	}

	//Set elements of A(2m, 3) matrix
	if (me_scheme == ScaleShiftsScheme)
	{
		A.replace(E, 0, 1);
		A.replace(E, m1, 2);
	}

	//Set elements of A(2m, 4) matrix
	else if (me_scheme == SimilarityScheme)
	{
		A.replace(-1.0 * Q(0, m1 - 1, 1, 1), 0, 1);
		A.replace(Q(0, m1 - 1, 0, 0), m1, 1);

		A.replace(E, 0, 2);
		A.replace(E, m1, 3);
	}

	//A.print();

	//Set initial weights W = I
	W = eye(2 * m1, 2 * m1, 1.0);

	//Create Y vector
	Matrix <T> Y(2 * m1, 1);
	Y.replace(P(0, m1 - 1, 0, 0), 0, 0);
	Y.replace(P(0, m1 - 1, 1, 1), m1, 0);

	//Compute initial residuals
	Matrix <T> beta = pinv1(trans(A) * W * A) * trans(A) * W  * Y;

	//Compute initial residuals
	Eps = Y - A * beta;
	f_init = norm(trans(Eps) * W * Eps);

	//Set parameters
	T df = 2 * tol;
	iter = 0;
	I = ones(2 * m1, 1, 1);

	//Main IRLS procedure
	while ((df > tol) && (iter < max_iter))
	{
		//Compute new beta using the least squares: more accurate than a common inversion
		beta = pinv1(trans(A) * W * A) * trans(A) * W  * Y;

		//Compute residuals
		Eps = Y - A * beta;

		//Compute median absolute deviation
		const T Sn = mad(Eps, 0);

		/*/
		//Cofactor matrix
		Matrix <T> Ql = inv(W);

		//Covariance matrix
		Matrix <T> Qv = diag(diag(inv(W) - A * inv(trans(A) * W * A) * trans(A)));

		//Covariance matrix
		Matrix <T> Cl = Sn * Sn* Ql;
		*/

		//Apply weight functions
		for (unsigned int i = 0; i < 2 * m1; i++)
		{
			//Compute normalized residual
			const T epsi = abs(Eps(i, 0) / Sn);

			//Huber estimator
			if (me_weight_function == HuberFunction)
			{
				if (epsi <= k)
				{
					W(i, i) = 1;
				}

				else
				{
					W(i, i) = k / epsi;
				}
			}

			//Andrew function
			else if (me_weight_function == AndrewFunction)
			{
				if (epsi <= k * M_PI)
				{
					W(i, i) = k / epsi * sin( epsi / k );
				}

				else
				{
					W(i, i) = 0;
				}
			}

			//Tukey biweight function
			else if (me_weight_function == TukeyFunction)
			{
				if (epsi <= k)
				{
					const T epsi2 = epsi * epsi / (k * k);
					W(i, i) = (1 - epsi2) * (1 - epsi2);
				}

				else
				{
					W(i, i) = 0;
				}
			}

			//Danish method (WLS28)
			else if (me_weight_function == DanishFunction)
			{
				if (epsi < k)
				{
					W(i, i) = 1;
				}

				else
				{
					W(i, i) = exp(-epsi * epsi / k);
				}
			}

			//Modified Danish method(WLS23, WLS15, WLS26)
			else if (me_weight_function == DanishFunction2)
			{
				if (epsi > k)
				{
					W(i, i) *=  exp( -epsi / k );
				}
			}

			//IGCIII, Yang(WLS23, WLS24)
			else if (me_weight_function == YangFunction)
			{
				const T c0 = 1.345;
				const T c1 = 3.0;

				//Good
				if (epsi <= c0)
				{
					W(i, i) = 1;
				}

				//Suspected
				else if ((epsi > c0) && (epsi <= c1))
				{
					const T epsin = (c1 - epsi) / (c1 - c0);
					W(i, i) = c0 / epsi * epsin * epsin;
				}

				//Outliers
				else if (epsi > c1)
				{
					W(i, i) = 0;
				}
			}
		}

		//Compute residuals
		const T f2 = norm(trans(Eps) * W * Eps);

		//Compute the condition
		df = fabs(f2 - f);

		//Assign as old value
		f = f2;

		//Increment itrerations
		iter++;

		const Matrix <T> WD = diag(W);
		//WD.print();
	}

	//Find indices of possible outliers
	unsigned int index = 0;
	for (unsigned int i = 0; i < 2 * m1; i++)
	{
		//Observation = outlier ?
		if (((me_weight_function == HuberFunction) && (W(i, i) < 1)) || ((me_weight_function == AndrewFunction) && (W(i, i) < wmax)) 
		 || ((me_weight_function != HuberFunction) && (me_weight_function != AndrewFunction) && (W(i, i) <= wmin)))
		{
			I(index, 0) = 0;
		}

		index = index + 1;
	}

	//I.print();
	//Matrix <T> WW = diag(W);
	//WW.print();
}


template <typename Point1, typename Point2, typename TKey>
void Outliers::findOutliersIRLS(const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source,
	Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest, TKey & min_key, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs)
{
	//Removing outliers using the iteratively IRLS
	//Outliers are dectected using weighted transformation and Danish method
	//Method published by Wieser and Brunner, 2000
	const unsigned int n_global_points_source = global_points_source.size();
	typename Point1::Type MIN_DEV_DIFFERENCE = 0.0001, std_dev_old = 0;

	//Accuracy parameters
	TAccuracyCharacteristics <typename Point1::Type> accuracy;

	//Throw exception: not enough points
	if (n_global_points_source < 3)
	{
		throw BadDataException("BadDataException: not enough global points to be optimized. ", "Can not find optimal transformation key.");
	}

	//Initialize weights
	typename TWeights <typename Point1::Type> ::Type weights(n_global_points_source, 1);

	//First iteration: compute as unweighted triangulation
	Container <Point1 *> transformed_points_init;
	Transformation2D::getTransformKey(global_points_source, local_points_source, weights, min_key);
	Transformation2D::transform(global_points_source, local_points_source, transformed_points_init, min_key);
	accuracy = Transformation2D::getAccuracyCharacteristics(global_points_source, local_points_source, transformed_points_init, min_key, weights);

	//Compute initial q_ee
	typename TItemsList <typename Point1::Type>::Type q_ee(n_global_points_source, 1), q_ee_inv(n_global_points_source, 1);

	for (unsigned int i = 0; i < n_global_points_source; i++)
	{
		q_ee[i] = (1.0 / weights[i] - accuracy.q_xx[i]);
	}

	//Initialize accuracy
	accuracy.std_dev = MAX_FLOAT;

	//Perform iterations until no change of standard deviation
	do
	{
		//Remember old deviation
		std_dev_old = accuracy.std_dev;

		//Invert Q_ee matrix
		for (unsigned int i = 0; i < n_global_points_source; i++)
			q_ee_inv[i] = 1.0 / q_ee[i];

		//Compute weighted transformation
		Container <Point1 *> transformed_points;
		Transformation2D::getTransformKey(global_points_source, local_points_source, q_ee_inv, min_key);
		Transformation2D::transform(global_points_source, local_points_source, transformed_points, min_key);

		//Compute accuracy characteristics for weighted transformation
		accuracy = Transformation2D::getAccuracyCharacteristics(global_points_source, local_points_source, transformed_points, min_key, q_ee_inv);

		//Remove outliers using Danish method
		for (unsigned int i = 0; i < n_global_points_source; i++)
		{
			//Compute new Q_ee
			q_ee[i] = (1.0 / weights[i] - accuracy.q_xx[i]);

			//Danish method: decrease weight for outliers
			if ((accuracy.res[i].res_xy > 2.0 * accuracy.std_dev * sqrt(q_ee[i])) /*&& (accuracy.res[i].res_xy > 200 )*/)
			{
				weights[i] *= exp(-(accuracy.res[i].res_xy / (2.0 * accuracy.std_dev * sqrt(q_ee[i]))));
			}
		}

	} while (fabs(std_dev_old - accuracy.std_dev) > MIN_DEV_DIFFERENCE);

	//Create best pairs and add them to the key: removing outliers having decreased weights
	for (unsigned int i = 0; i < n_global_points_source; i++)
	{
		if (weights[i] == 1.0)
		{
			min_pairs.push_back(std::make_pair(accuracy.res[i].res_xy, i));
		}
	}

	//Rearrange points: add best points to both dest lists
	Transformation2D::rearrangePoints(global_points_source, local_points_source, global_points_dest, local_points_dest, min_pairs);

	//Compute transformation key using non-weighted key
	Transformation2D::getTransformKey(global_points_dest, local_points_dest, min_key);
}


template <typename Point1, typename Point2>
void Outliers::findOutliersME(const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source, Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest, const typename Point1::Type k, const typename Point1::Type tol,
	const TMEstimatorsScheme me_scheme, TMEstimatorsWeightFunction me_weight_function, const unsigned int max_iter, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, typename Point1::Type &f_init, typename Point1::Type &f, unsigned int &iter)
{
	//Removing outliers using the M-estimated
	unsigned int m = global_points_source.size();

	//Convert points to matrices
	Matrix <typename Point1::Type> P(m, 2), Q(m, 2), W = eye(2 * m, 2 * m, 1.0), Eps( 2 *m , 1);
	Matrix <unsigned int> I(2 * m, 1);

	for (unsigned int i = 0; i < m; i++)
	{
		P(i, 0) = global_points_source[i]->getX(); P(i, 1) = global_points_source[i]->getY();
		Q(i, 0) = local_points_source[i]->getX(); Q(i, 1) = local_points_source[i]->getY();
	}
	
	//P.print();
	//Q.print();

	//Find outliers using M-estimators
	f_init = 0, f = 0;
	findOutliersME(P, Q, k, tol, me_scheme, me_weight_function, max_iter, W, I,Eps, f_init, f, iter);

	Matrix <double> WD = diag(W);
	WD.print(); 

	//Create pair: residual, point number, add correct points to the list
	for (unsigned int i = 0; i < m; i++)
	{
		//Observation is not an outlier
		if ((I(i, 0) == 1.0) && (I(i + m, 0) == 1.0))
		{
			//Compute residual
			const typename Point1::Type res_xy = sqrt(Eps(i, 0) * Eps(i, 0) + Eps(i + m, 0) * Eps(i + m, 0));

			//Add to the correct pairs
			min_pairs.push_back(std::make_pair(res_xy, i));

			//Create copies of correct points
			global_points_dest.push_back(new Point1(global_points_source[i]));
			local_points_dest.push_back(new Point2(local_points_source[i]));
		}
	}

	I.print();
}


template <typename T>
void Outliers::createKBestPairsOfPoints(const TAccuracyCharacteristics <T> &deviations, typename TDevIndexPairs <T>::Type & point_pairs, const float perc_ratio)
{
	//Create pairs of all points, sort by standard deviation and erase n-k points with greatest outliers
	const unsigned int n = deviations.res.size();

	//Clear old pairs
	point_pairs.clear();

	//Create new pairs
	for (unsigned int j = 0; j < n; j++)
	{
		point_pairs.push_back(std::make_pair(deviations.res[j].res_xy, j));
	}

	//Sort pairs by standard deviation
	std::sort(point_pairs.begin(), point_pairs.end(), sortPointPairsByResiduals <T>());

	//Let only first k-items
	point_pairs.erase(point_pairs.begin() + n * perc_ratio, point_pairs.end());
}


#endif
