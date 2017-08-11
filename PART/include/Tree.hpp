#ifndef H_TREE_HPP
#define H_TREE_HPP

#include <armadillo>
#include <vector>

using namespace arma;

class RawMCMC;
class MCestimate;

/*! @file Tree.cpp
@brief Builds a binary tree in which each leaf node represents an element of the partition of the input sampling area (stored in input RawMCMC object).

The constructor of this class calls the private method Tree::build() that recursively determines a partition of the current sampling area.
Each element of \c leaves vector represents a partition's block and stores information on the corresponding MCMC samples.
Since only leaf nodes are relevant for posterior resampling, outside this class no information is available on those nodes in the binary tree that do not represent leaf nodes.
This class is used in oneStageResampling().

@see RawMCMC for additional information on how to handle MCMC samples and sampling area information.
@see MCestimate for additional information on relevant statistics for MCMC samples.
@see Tree::build() for additional information on how to grow the binary tree.
@see Tree::Node() for additional information on how each node in the bunary tree is represented.
@see oneStageResampling() for additional information on how the values stored in \a leaves vector will be used for posterior resampling.
*/
class Tree{

	public:

		//! Constructor of binary tree (represented through \c leaves vector)
		/*!
		@param kdCut (\c uword): \c true if cutting point determined as median of MCMC samples (KD-cut); \c false if cutting point determined as maximizer of empirical log-likelihood (ML-cut).
		@param minNumPointsBlock (\c uword): minimum number of MCMC samples falling into each partition's block.
		@param minCutLength (\c double): minimum length of each dimension of partition's block.				
		*/
		Tree(const bool kdCut, const uword minNumPointsBlock, const double minCutLength);

		//! Method to define how to grow binary tree
		/*!
		@param rawMCMC (\c RawMCMC): RawMCMC object storing information on MCMC samples, corresponding subsets' indices and sampling area.
		@param dim_cut (\c uvec ): dimensions of sampling area that do not violate \a minCutLength.
		@param verbose (\c bool): \c true if all messages must be displayed; \c false if only basic information displayed.
		@param leaves (\c vector \c < \c MCestimate \c >): information regarding binary tree; each leaf node in \a leaves represents one block of the partition of the sampling area.
		
		This method recursively calls the private method build() to determine a partition of current sampling area.
		When the tree has been grown, the information on normalized and unnormalized densities, stored in each leaf node in \c leaves, are updated according to normalize().

		@see RawMCMC for additional details on how to handle information regarding MCMC samples.		
		@see MCestimate for additional details on how to store statistics regarding MCMC samples.
		@see Tree::build() for additional details on how to grow the partition tree.
		@see Tree::normalize() for additional details on how to normalize the estimated densities.
		*/
		void grow(const RawMCMC & rawMCMC, const uvec & dim_cut, const bool verbose, std::vector<MCestimate> & leaves);
		
	private:
		
		bool kd_cut; /*!< Defines method to be used to find cutting point; if \c true use median (KD), else use \c ml_cutting_rule() (ML). */
		uword min_num_points_block; /*!< Minimum number of samples falling into each leaf node. */
		double min_cut_length; /*!< Minimum length of sampling area along each dimension. */

		//! getter of \c kd_cut attribute
		const bool get_kd_cut() const;
		
		//! getter of \c min_num_points_block attribute
		const uword get_min_num_points_block() const;
		
		//! getter of \c min_cut_length attribute
		const double get_min_cut_length() const;

		//! setter of \c kd_cut attribute
		void set_kd_cut(const bool kdCut);
		
		//! setter of \c min_num_points_block attribute
		void set_min_num_points_block(const uword minNumPointsBlock);
		
		//! setter of \c min_cut_length attribute
		void set_min_cut_length(const double minCutLength);
		
		//! Struct to represent a generic node of partition tree.
		/*!
		Each \c Node object stores information regarding:
		- its relatives (pointers to parent, left and right child),
		- indices of samples falling into current node (w.r.t. samples matrix of parent node)
		- sampling area
		- dimensions along which search of valid cutting point could continue.
		\n This class is fundamental to grow the binary tree and is used in Tree::build().
		
		@see Tree::build() for additional information on how to grow the partition tree.
		*/
		struct Node{
			
			Node * parent; /*!< Pointer to parent node. */
			Node * childLeft; /*!< Pointer to left-child node. */
			Node * childRight; /*!< Pointer to right-child node. */
			uvec indices; /*!< Indices of samples falling into current node. */
			mat area; /*!< Sampling area. */
			uvec dim_cut; /*!< Dimensions of sampling area that do not violate min_cut_length. */

			//! Constructor of root node
			/*! Initialization of parent, childLeft, childRight attributes to nullptr */
			Node();
			
			//! Method to define attibutes of a child node and to link it to current node
			/*!
			@param dim (\c uword): dimension along which valid cutting point has been found.
			@param cut (\c double): value of cutting point.
			@param min_cut_length (\c double): minimum length of each dimension of partition's block.					
			@param is_left (\c bool): \c true ( \c false ) if input node need to be associated to current node as left (right) child.
			@param idx (\c uvec): subset of indices of samples (w.r.t. current node) that need to be stored in child node.
			@param node (\c Node): node to be linked as child of current node.
			*/	
			void generateChild(const uword dim, const double cut, const double min_cut_length, const bool is_left, const uvec & idx, Node & node);
			
		};
		
		//! Method to search for partition of current block.
		/*!
		@param rawMCMC (\c RawMCMC): RawMCMC object storing information on MCMC samples.
		@param node (\c Node): current node to be analyzed.
		@param cut (\c double): value of cutting point.
		@param dim (\c uword): dimension along which valid cutting point has been found.
		@param samples_dim (\c vec): column of matrix representing projection of samples belonging to current node along selected dimension \c dim.
		@param valid_cut (\c bool): \c true (\c false ) if valid cutting point has (not) been found.
		
		To search for a valid partition, the following steps are repeated:\n
		- randomly select one of the dimensions of the node's area along which a valid cut should be searched,
		- determine the candidate cutting point (according to \c kd_cut attribute) along the previously selected dimension,
		- call checkCut() to check the validity of previously determined candidate cutting point (i.e. for each resulting partition's block: area of block not too small and not too few samples falling into each block).
		
		The search stops when one of the following mutually exclusive conditions is met:
		- a valid cutting point has been found
		- all dimensions of the sampling area have been analyzed without finding any valid partition.
		
		@see RawMCMC for additional details on how to handle information regarding MCMC samples.
		@see Node for additional details on how to represent nodes in tree.
		@see ml_cutting_rule() for additional details on how to determine cutting point if \c kd_cut is \c false.
		@see checkCut() for additional details on how to establish the validity of a candidate cutting point.
		*/
		void findCut(const RawMCMC & rawMCMC, const Node & node, double & cut, uword & dim, vec & samples_dim, bool & valid_cut);
		
		//! Method to establish if candidate cutting point cut is valid
		/*!
		@param cut (\c double): value of cutting point.
		@param samples_dim (\c vec): column of samples belonging to current node and corresponding to selected dimension \c dim.
		@param l (\c double): left extreme of current sampling interval.
		@param r (\c double): right extreme of current sampling interval.
		
		@return valid_cut (\c bool): \c true (\c false ) if valid cutting point has (not) been found. 
		
		The candidate cutting point \c cut is valid only if both the following conditions hold:\n
		- area of each induced partition's block is not too small (i.e. does not violate \c min_cut_length)
		- number of samples falling into each block is not too small (i.e. does not violate \c min_num_points_block).
		
		@see findCut() for additional details on how to search for a candidate cutting point.
		*/
		bool checkCut(const double cut, const vec & samples_dim, const double l, const double r);
		
		//! Method to determine cutting point along a fixed dimension using Maximum Likelihood criterion
		/*!
		@param rawMCMC (\c RawMCMC): RawMCMC object storing information on MCMC samples.
		@param node (\c Node): current node to be analyzed.
		@param dim (\c uword): dimension along which cutting point must be searched.
		
		@return cut (\c double): value of the cutting point that could be used to partition the interval \c [ \c area \c(0,\c dim \c), \c area \c(1,\c dim \c) \c).
		
		According to this criterion, the best cutting point \c cut that defines a partition of current interval \c [ \c area \c(0,\c dim \c), \c area \c(1,\c dim \c) \c) is the sample along dimension \c dim such that the following objective function is maximized:
		\f[ \sum_{m=0}^{M-1} ( nLeft(m)*( \log(nLeft(m)) - \log(cut-area(0,dim))) ) + \sum_{m=0}^{M-1} ( nRight(m)*( \log(nRight(m)) - \log(area(1,dim)-cut)) ) \f]
		\n where:\n
		- nLeft(m) = number of MCMC samples in current node along dimension \c dim that belong to subset \a m and fall into \c [ \c area \c(0,\c dim \c), \c cut)
		- nRight(m) = number of MCMC samples in current node along dimension \c dim that belong to subset \a m and fall into \c [ \c cut, \c area \c(1,\c dim \c) \c)
		
		@see RawMCMC for additional details on how to handle information regarding MCMC samples.
		@see findCut() for additional details on how to search for a candidate cutting point.
		*/		
		double ml_cutting_rule(const RawMCMC & rawMCMC, const Node & node, const uword dim);
		
		//! Method to recursively build binary tree
		/*!
		@param rawMCMC (\c RawMCMC): RawMCMC object storing information on MCMC samples.
		@param node (\c Node): current node to be analyzed.
		@param verbose (\c bool) : \c true if all messages must be displayed; \c false if only basic information displayed.
		@param leaves (\c vector \c < \c MCestimate \c >): information regarding binary tree; each leaf node in \c leaves represents one partition's clock of sampling area.
		
		For given input \c node, the function searches for a valid partition of the current sampling area calling method findCut().	
		
		In case a valid cut has been found, the current node does not represent yet an element of the sampling area's partition. \n
		In such a case, the search should proceed calling build() onto children nodes, each of them generated through Node::generateChild.

		In case all dimensions of the sampling area have been analyzed and no valid cutting point has been found, the current node represents an element of the sampling area's partition.\n
		Since we will need to use the information on samples falling in current node (represented through an object of class MCestimate), we store them in \c leaves vector.
		
		The result of the overall execution of the method build() is to grow the \c leaves vector, storing information regarding all those blocks that represent a partition of the input sampling area.
		
		@see findCut() for additional details on characteristics of valid partition.
		@see MCestimate for additional details on how to store statistics regarding MCMC samples.
		*/
		void build(const RawMCMC & rawMCMC, Node & node, const bool verbose, std::vector<MCestimate> & leaves);
		
		//! Method to normalize the binary tree (represented through \c leaves vector)
		/*!
		@param leaves (\c vector \c < \c MCestimate \c >): information regarding current tree, in which each leaf node represents one element of the partition of the sampling area.
		
		This method is used only if MCMC samples have been drawn from more than one subset.
		
		Indeed, if M=1, the estimate of the posterior density is given by the K-block histogram corresponding to subset 0 hence:
					\f[ dens(\theta|data) = \sum_{k=0}^{K-1} ( n_{block}(0,k) / N(0) * 1 / |block(k)| * g_k(\theta) ) \f]
			  where\n
					- K = leaves.size(): number of leaves in current partition tree\n
					- n_block(0,k): number of samples drawn from subset m=0 falling into k-th leaf node\n
					- N(0): number of samples drawn from subset m=0\n
					- block(k): partition of sampling area at k-th leaf node\n
					- g_k(\theta): local density based on samples' estimates at k-th leaf node\n
						i.e. either g_k(\theta) = Unif(block(k)) where block(k) is partition of sampling area at k-th leaf node or g_k(\theta) = Gauss(mean(k), cov(k)) where mean(k), cov(k) are mean and covariance estimates deriving from k-th leaf node.\n
		
		If the number of subsets is more than one, the weight of each mixture component need to be normalized by global normalization constant.
		The normalization constant is determined as the sum of the unnormalized weight at each partition's block.\n
		If M>1, the estimate of the posterior density need to be normalized as follows:
					\f[ dens(\theta|data) = \sum_{k=0}^{K-1} ( prob(k) * g_k(\theta) ) \f]
			  where\n
					- K = leaves.size(): number of leaves in current partition tree\n
					- prob(k) = prob_hat(k) / Z: weight of k-th mixture component\n
					where\n
						\f[ \hat{prob}(k) = (1 / |block(k)|)^(M-1) * prod_{m=0}^{M-1} ( n_{block}(m,k) / N(m) ) \f]
						is the unnormalized weight of k-th weight component, with: \n
							- n_block(m,k): number of samples drawn from subset m falling into k-th leaf node\n
							- N(m): number of samples drawn from subset m\n
							- block(k): partition of sampling area at k-th leaf node\n
						Z = sum_{k=0:K-1} (prob_hat(k)): normalizing constant\n
					- g_k(\theta): local density based on samples' estimates at k-th leaf node\n
						i.e. either g_k(\theta) = Unif(block(k)) where block(k) is partition of area at k-th leaf node or g_k(\theta) = Gauss(mean(k), cov(k)) where mean(k), cov(k) are mean and covariance estimates deriving from k-th leaf node
		
		@see MCestimate for additional details on how to store statistics regarding MCMC samples.
		*/
		void normalize(std::vector<MCestimate> & leaves);
};
#endif