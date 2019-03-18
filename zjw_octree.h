#pragma once

/*!
 * \file zjw_octree.h
 *
 * \author King
 * \date 四月 2018
 *  this head file and the implements c++ souce code has set the basis octree build and opration
 *  interface. which implements in real application should look at zjw_pcs_octree.h
 */

#include <assert.h>
#include "util/zjw_macro.h"
#include <util/zjw_math.h>


/*
* 输入：
*	bit: 表示对x,y,z中的哪一维度进行判断
*   p: 表示的当前需要进行判断的点云帧中的顶点
*   mid: 表示的是当前Node的中点的坐标，需要根据这个坐标分出8个象限
* 输出：
*	i: 表示当前p顶点，应该在当前的node上的哪一个象限
*   newMin: 表示新的象限边界的最小点的值
*   newMax: 表示新的象限边界的最大点的值
*/
#define COMPUTE_SIDE(i, bit, p, mid, newMin, newMax) \
if (p >= mid)         \
{                     \
    i |= bit;         \
    newMin = mid;     \
}                     \
else                  \
{                     \
    newMax = mid;     \
}


 //定义了节点数据类型为 N 的八叉树
template <class N>
class Octree
{
public:
	//------------define 八叉树的结点 use in the tree. ---------------
	class OctreeNode
	{
	public:
		//这个节点的数据
		N nodeData;
		//非叶子节点是-1，其他的从0,1,2,3一次编号
		int leafFlag;
		//min的顶点和Max的顶点确定 节点的空间范围
		Vec3 min;
		Vec3 max;
		OctreeNode* children[8];

	public:
		OctreeNode()
		{
			leafFlag = -1;
			//初始化八个孩子的结点
			for (int i = 0; i < 8; i++)
				children[i] = 0;
		}

		virtual ~OctreeNode()
		{
			for (int i = 0; i < 8; i++)
				if (children[i])
					delete children[i];
		}
	};

	//min的顶点和Max的顶点确定 节点的空间范围
	Vec3 min;
	Vec3 max;
	// cellSize is the minimum subdivided cell size.
	Vec3 cellSize;
	OctreeNode* root;

public:
	Octree(Vec3 min, Vec3 max, Vec3 cellSize) : min(min), max(max), cellSize(cellSize){
		root = new OctreeNode();
		//min的顶点和Max的顶点确定 节点的空间范围
		root->min = min;
		root->min = max;			
	}
	virtual ~Octree() { delete root; }

	//遍历所有的Octree中做操作的一个类。
	class Callback
	{
	public:
		// Return value: true = continue; false = abort.
		virtual bool operator()(const Vec3 min, const Vec3 max, N& nodeData) { return true; }

		// Return value: true = continue; false = abort.
		virtual bool operator()(const Vec3 min, const Vec3 max, OctreeNode* currNode) { return true; }
	};

	//根据深度遍历的方法，把点pos插入到八叉树的叶子节点中，直到叶子节点的size小于 cellSize.
	//N& getCell(const float pos[3], Callback* callback = NULL)
	N& getCell(const Vec3 ppos, Callback* callback = NULL)
	{
		//Point3f ppos(pos);
#ifdef ZJW_DEBUG
//		cout << (ppos >= min) << endl;
	//	cout << (ppos < max) << endl;
#endif
		assert(ppos >= min && ppos < max);
		Vec3 currMin(min);
		Vec3 currMax(max);
		Vec3 delta = max - min;
		if (!root)
			root = new OctreeNode();
		OctreeNode* currNode = root;
		while (delta >= cellSize)
		{
			bool shouldContinue = true;
			if (callback)
				shouldContinue = callback->operator()(currMin, currMax, currNode->nodeData);
			if (!shouldContinue)
				break;
			Vec3 mid = (delta * 0.5f) + currMin;
			Vec3 newMin(currMin);
			Vec3 newMax(currMax);
			int index = 0;
			//判断ppos这个点，在当前的这个node范围中，应该属于哪个象限。
			COMPUTE_SIDE(index, 1, ppos.x, mid.x, newMin.x, newMax.x)
			COMPUTE_SIDE(index, 2, ppos.y, mid.y, newMin.y, newMax.y)
			COMPUTE_SIDE(index, 4, ppos.z, mid.z, newMin.z, newMax.z)

			if (!(currNode->children[index]))
				currNode->children[index] = new OctreeNode();

			currNode = currNode->children[index];
			//每个节点保留自己的范围
			currNode->min = newMin;
			currNode->max = newMax;

			currMin = newMin;
			currMax = newMax;
			delta = currMax - currMin;
		}
		return currNode->nodeData;
	}

	//递归遍历整个树.并进行某种操作。
	void traverse(Callback* callback)
	{
		assert(callback);
		traverseRecursive(callback, min, max, root);
	}

	//删除整个树，释放所有节点的空间
	void clear()
	{
		delete root;
		root = NULL;
	}

	//遍历的内部类
	class Iterator
	{
	public:
		Iterator getChild(int i)
		{
			return Iterator(currNode->children[i]);
		}
		N* getData()
		{
			if (currNode)
				return &currNode->_nodeData;
			else return NULL;
		}
	protected:
		//当前node的指针
		OctreeNode* currNode;

		Iterator(OctreeNode* node) : currNode(node) {}
		friend class Octree;
	};

	Iterator getIterator()
	{
		return Iterator(root);
	}

protected:

	//递归遍历整个树
	void traverseRecursive(Callback* callback, const Vec3& currMin, const Vec3& currMax, OctreeNode* currNode)
	{
		//为空返回
		if (!currNode)
			return;
		//判断是否符合其他终止条件
		bool shouldContinue = callback->operator()(currMin, currMax, currNode);
		if (!shouldContinue)
			return;

		//mid 得到cell中间的点
		Vec3 delta = currMax - currMin;
		Vec3 mid = (delta * 0.5f) + currMin;
		//分成子八个节点
		traverseRecursive(callback, currMin, mid, currNode->children[0]);
		traverseRecursive(callback, Vec3(mid.x, currMin.y, currMin.z), Vec3(currMax.x, mid.y, mid.z), currNode->children[1]);
		traverseRecursive(callback, Vec3(currMin.x, mid.y, currMin.z), Vec3(mid.x, currMax.y, mid.z), currNode->children[2]);
		traverseRecursive(callback, Vec3(mid.x, mid.y, currMin.z), Vec3(currMax.x, currMax.y, mid.z), currNode->children[3]);
		traverseRecursive(callback, Vec3(currMin.x, currMin.y, mid.z), Vec3(mid.x, mid.y, currMax.z), currNode->children[4]);
		traverseRecursive(callback, Vec3(mid.x, currMin.y, mid.z), Vec3(currMax.x, mid.y, currMax.z), currNode->children[5]);
		traverseRecursive(callback, Vec3(currMin.x, mid.y, mid.z), Vec3(mid.x, currMax.y, currMax.z), currNode->children[6]);
		traverseRecursive(callback, mid, currMax, currNode->children[7]);
	}
};
