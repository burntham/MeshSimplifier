/*
 * simplifier.h
 *
 *  Created on: Sep 17, 2014
 *      Author: brunt
 *  This is an abstract class used to structure the simplifiers
 */

#include <iostream>

#ifndef SIMPLIFIER_H_
#define SIMPLIFIER_H_

namespace brndan022 {

template <typename M> class Simplifier {
public:
	virtual void simplify(M&) =0;
	virtual void setParameters(int,char **) =0;
	virtual ~Simplifier(){};
};

} /* namespace brndan022 */

#endif /* SIMPLIFIER_H_ */
