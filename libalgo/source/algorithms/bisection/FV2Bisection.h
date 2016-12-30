
#ifndef FV2BISECTION_H
#define FV2BISECTION_H

template <typename T>
//template <typename T, typename FunctionV>
class FV2Bisection
{
	private:
		Matrix <T> &X;

	public:
		void operator() ( const T lon0, T &cost)
		{
			X(4, 0) = lon0;
			Matrix <T> Y(1,1), W(1,1), V(1,1);
			//FunctionV(X, Y, V, W);	
			FAnalyzeProjV2(X, Y, V, W);

			cost = norm(trans(W)*W);
		}	
	
};

#endif