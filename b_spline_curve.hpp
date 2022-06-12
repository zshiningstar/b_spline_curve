#pragma once

#include<iostream>
#include<vector>
#include<cmath>
#include<string>

/* @brief 点的例子
 */
class Point
{
public:
	double x,y,z;
	Point(){}
	Point(double _x, double _y, double _z=0.0):
		x(_x), y(_y), z(_z){}
	
	float disTo(const Point& point) const 
	{
		float dx = point.x - x;
		float dy = point.y - y;
		return sqrt(dx*dx + dy*dy);
	}

    float disTo(Point* const point) const
    {
        float dx = point->x - x;
        float dy = point->y - y;
        return sqrt(dx*dx + dy*dy);
    }
};

template <class TYPE>
void showVector(std::vector<TYPE>& vec) {
	for (int i = 0; i < vec.size(); i++) {
		std::cout << vec[i] << " ";
	}
	std::cout << std::endl;
}

/* @brief linspace函数
 */
std::vector<double> linSpace(double x1, double x2, int n){
    std::vector<double> ans;
    int i = 0;
    double d = (x2-x1) / (n-1);
    for(i = 0; i < n; i++){
        double temp_point = x1 + i * d;
        ans.emplace_back(temp_point);
    }
    
    return ans;
}
/* @brief u_quasi_uniform函数
 */
std::vector<double> u_quasi_uniform(int n, int k){
    //准均匀B样条的节点向量计算，共n+1个控制顶点，k次B样条，k+1阶
    std::vector<double> ans(n+k+2,0);
    //曲线段数
    int piecewise = n - k + 1;
//    std::cout << "piecewise:    " << piecewise << std::endl;
    if(piecewise == 1){
        for(int i = k+2; i <= n+k+2; i++){
            ans[i] = 1.0;
        }
//        std::cout << "piecewise:   " << piecewise << std::endl;
    }else{
        int flag = 0;
        while(flag != piecewise){
            ans[k+flag+1] = ans[k+flag] + 1.0/piecewise;
            flag += 1;
            showVector(ans);
        }
        for(int i = n+2; i <= n+k+2; i++){
            ans[i] = 1.0;
        }
    }
    return ans;
}

double baseFunction(int i, int k, double u, std::vector<double> nodeVector){
    double bik = 0.0;
    if(k == 0){
        if((u >= nodeVector[i]) && (u < nodeVector[i+1])){
            bik = 1.0;
        }else{
            bik = 0.0;
        }
    }else{
        double length1 = nodeVector[i+k] - nodeVector[i];
        double length2 = nodeVector[i+k+1] - nodeVector[i+1];
//        std::cout << "length1:    " << length1 << "   " << "length2   " << length2 << std::endl;
        if(length1 == 0){//规定0/0=0
            length1 = 1;
        }
        if(length2 == 0){
            length2 = 1;
        }
        bik = (u - nodeVector[i]) / length1 * baseFunction(i, k-1, u, nodeVector) + 
              (nodeVector[i+k+1] - u) / length2 * baseFunction(i+1, k-1, u, nodeVector);
    }
    return bik;
}


/* @brief B样条曲线插值
 * @P:控制点,6个控制点满足曲率连续
 * @tpye:B样条类型
 * @k:曲线为k阶,k-1次
 * @path:返回值
 */
std::vector<Point> generateBsplinePath(std::vector<Point>& P, const bool type){
    int k = 4;
    int n = P.size() - 1;
    std::vector<Point> path;
    std::vector<double> bik;
    if(type){//均匀B样条
        std::vector<double> nodeVector = linSpace(0, 1, n + k + 1);
        showVector(nodeVector);
        double u = (k-1.0) / (n+k+1.0);
        for(u = (k-1.0) / (n+k+1.0); u <= (n+2.0) / (n+k+1.0) + 0.001; u += 0.001){
            bik.clear();
            for(int i = 0; i <= n; i++){
                double temp_bik = baseFunction(i, k-1, u, nodeVector);
                bik.emplace_back(temp_bik);
            }
            showVector(bik);
            Point pu(0,0);
            for(int i = 0; i <= n; i++){
                pu.x = pu.x + P[i].x * bik[i];
                pu.y = pu.y + P[i].y * bik[i];
            }
            path.emplace_back(pu);
        }
    }else{//准均匀B样条
        std::vector<double> nodeVector = u_quasi_uniform(n, k-1);
        showVector(nodeVector);
        double u = 0.0;
        for(; u <= 1-0.05; u += 0.005){
            bik.clear();
            for(int i = 0; i <= n; i++){
                double temp_bik = baseFunction(i, k-1, u, nodeVector);
                bik.emplace_back(temp_bik);
            }
//            showVector(bik);
            Point pu(0,0);
            for(int i = 0; i <= n; i++){
                pu.x = pu.x + P[i].x * bik[i];
                pu.y = pu.y + P[i].y * bik[i];
            }
            path.emplace_back(pu);
        }
    }
    return path;
}
