% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/method.R
\name{getArgumentDefaults}
\alias{getArgumentDefaults}
\alias{getArgumentDefaults,lcMethod-method}
\title{Default argument values for the given method specification}
\usage{
getArgumentDefaults(object, ...)

\S4method{getArgumentDefaults}{lcMethod}(object)
}
\arguments{
\item{object}{The method specification object.}

\item{...}{Not used.}
}
\value{
A \verb{named list} of argument values.
}
\description{
Returns the default arguments associated with the respective \code{lcMethod} subclass.
These arguments are automatically included into the \code{lcMethod} object during initialization.
}
\section{Implementation}{

Although implementing this method is optional, it prevents users from
having to specify all arguments every time they want to create a method specification.

In this example, most of the default arguments are defined as arguments of the function
\code{lcMethodExample}, which we can include in the list by calling \link{formals}. Copying the arguments from functions
is especially useful when your method implementation is based on an existing function.
\preformatted{
setMethod("getArgumentDefaults", "lcMethodExample", function(object) {
  list(
    formals(lcMethodExample),
    formals(funFEM::funFEM),
    extra = Value ~ 1,
    tol = 1e-4,
    callNextMethod()
  )
})
}

It is recommended to add \code{callNextMethod()} to the end of the list.
This enables inheriting the default arguments from superclasses.
}

\seealso{
\link{getArgumentExclusions}

\link{lcMethod}

Other lcMethod implementations: 
\code{\link{getArgumentExclusions}()},
\code{\link{lcMethod-class}},
\code{\link{lcMethodAkmedoids}},
\code{\link{lcMethodCrimCV}},
\code{\link{lcMethodDtwclust}},
\code{\link{lcMethodFeature}},
\code{\link{lcMethodFunFEM}},
\code{\link{lcMethodFunction}},
\code{\link{lcMethodGCKM}},
\code{\link{lcMethodKML}},
\code{\link{lcMethodLMKM}},
\code{\link{lcMethodLcmmGBTM}},
\code{\link{lcMethodLcmmGMM}},
\code{\link{lcMethodMclustLLPA}},
\code{\link{lcMethodMixAK_GLMM}},
\code{\link{lcMethodMixtoolsGMM}},
\code{\link{lcMethodMixtoolsNPRM}},
\code{\link{lcMethodRandom}},
\code{\link{lcMethodStratify}}
}
\concept{lcMethod implementations}
