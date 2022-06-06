import importlib
import functools
import numpy as np
import itertools
import inspect
from typing import Callable, List

#for documentation generation
__arg_description = dict(
    N='polymer chain length',
    sigma='polymer brush grafting density',
    chi='Flory-Huggins parameter polymer-solvent',
    z='distance from grafting surface',
    R='distance from to grafting surface to the restriction',
    pore_Radius='grafted pore radius',
    ph="particle height",
    pw="particle width",
    pc="particle center",
    chi_PC='Flory-Huggins parameter polymer-particle',
    expansion_coefs="polynomial coefficients of gamma - phi dependency",
    eta="topological parameter",
    g="dendron generations",
    q="dendron functionality",
    percentile = "percentile to calculate D for smoothed phi profile"
)

#for documentation generation
__return_description = dict(
    D="polymer brush's thickness",
    osmotic_free_energy="osmotic part of the insertion free energy change",
    surface_free_energy="surface interaction contribution of the insertion free energy change",
    total_free_energy="total free energy of an inserted particle",
    phi='local polymer volume fraction',
    phi_D='local polymer volume fraction at the brush\'s end',
    Pi='local osmotic pressure',
    critical_pore_radius='minimal radius of an open pore'
)


def combine_signatures(
        callables: List[Callable],
        exclude=[],
        return_annotation=inspect.Signature.empty,
        positional_or_keyword_order : List[str] = None,
):

    signatures = [inspect.signature(c) for c in callables]
    params = {}
    default_params = {}
    for sig in signatures:
        for key, param in sig.parameters.items():
            if key in exclude:
                pass
            elif param.default == inspect.Parameter.empty:
                params[key] = param
                if key in default_params:
                    print(f"for '{key}' default value is inconsistent")
                    del default_params[key]
            else:
                default_params[key] = param
                if key in params:
                    print(f"for '{key}' default value is inconsistent")
                    del params[key]

    params.update(default_params)

    if positional_or_keyword_order:
        #if any(params[key].default != inspect.Parameter.empty for key in positional_or_keyword_order)
        ordered_params = {key : params[key] for key in positional_or_keyword_order if key in params}
        ordered_params.update(params)

        params = ordered_params
        
    
    signature = inspect.Signature(
        params.values(), return_annotation=return_annotation)

    return signature


def generate_docstring(signature: inspect.Signature, return_name):
    header = f'Calculates {__return_description[return_name]} for given arguments'

    args_desc = []
    for k, v in signature.parameters.items():
        if v.default is inspect.Parameter.empty:
            args_desc.append(
                f'{k} ({v.annotation.__name__}) : {__arg_description[k]}')
        else:
            args_desc.append(
                f'{k} ({v.annotation.__name__}, optional) : {__arg_description[k]}. Defaults to {v.default}')
    args_section = '\n\nArgs:\n' +\
        '\n'.join(args_desc)

    return_section = '\n\nReturns:' + \
        f'\n({signature.return_annotation.__name__}): ' + \
        f'{__return_description[return_name]}'

    return header+args_section+return_section


def generate_annotation(signature: inspect.Signature):
    annotations = {k: v.annotation for k, v in signature.parameters.items()}
    annotations["return"] = signature.return_annotation
    return annotations


def doc_decorator(
    callables : List[Callable],
    return_annotation=inspect.Signature.empty,
    return_name = None,
    exclude=[],
    #exclude_ignore_key_error = True,
    positional_or_keyword_order : List[str] = None
    ):
    signature = combine_signatures(
            callables,
            exclude, 
            return_annotation,
            positional_or_keyword_order
        )
    def decorator(func : Callable) -> Callable:
        def wrapped(*args, **kwargs):
            nonlocal signature
            #if any([kw not in signature.parameters for kw in kwargs]):
            #    raise KeyError("unexpected keyword argument")
            for kw in kwargs:
                if kw not in signature.parameters:
                    raise TypeError(f"'{kw}' unexpected keyword argument")
                if kw in exclude:
                    raise TypeError(f"'{kw}' keyword argument has been excluded")
            args_to_kwargs = {param : arg for param, arg in zip(signature.parameters, args)}
            return func(**args_to_kwargs, **kwargs) 
        wrapped.__signature__ = signature
        wrapped.__annotations__ = generate_annotation(signature)
        
        nonlocal return_name
        if return_name is None: return_name = func.__name__
        wrapped.__doc__ = generate_docstring(signature, return_name=return_name)
        wrapped.__name__ = func.__name__
        wrapped.__module__ = func.__module__
        return wrapped
    return decorator


def ignore_extra_kwargs(func : Callable) -> Callable:
    """Decorator that allows to call functions with redundant kwargs, 
    so that the function only uses required ones.
    """
    parameters = inspect.signature(func).parameters

    #if args_order is None:
    #    args_order = parameters.keys()
    #else:
    #    args_order = {k : parameters[k] for k in args_order if k in parameters}
    
    def wrapped(**kwargs):
        #new_args = {k : arg for k, arg in zip(args_order, args)}
        new_kwargs = {k: kwargs[k] for k in kwargs if k in parameters}
        return func(**new_kwargs)
    #wrapped.__annotations__ = func.__annotations__
    wrapped.__doc__ = func.__doc__
    wrapped.__name__ = func.__name__
    wrapped.__module__ = func.__module__
    return wrapped


def make_partial(*args, ignore_keyerror = False, **kwargs) -> Callable:
    def decorator(func):
        signature = inspect.signature(func)
        args_to_kwargs = {param : arg for param, arg in zip(signature.parameters, args)}
        all_args = dict(**args_to_kwargs, **kwargs)
       
        wrapped_params = {k : v for k, v in signature.parameters.items() if k not in all_args}

        def wrapped(*args2, **kwargs2):
            args_to_kwargs2 = {param : arg for param, arg in zip(wrapped_params, args2)}
            all_args2 = dict(**args_to_kwargs2, **kwargs2)
            if ignore_keyerror:
                return ignore_extra_kwargs(func)(**all_args, **all_args2)
            else:
                return func(**all_args, **all_args2)
        
        wrapped_signature = inspect.Signature(wrapped_params.values(), return_annotation=signature.return_annotation)
        #print(wrapped_signature)

        wrapped.__signature__ = wrapped_signature
        wrapped.__annotations__ = generate_annotation(wrapped_signature)
        wrapped.__doc__ = generate_docstring(signature, return_name=func.__name__)
        wrapped.__doc__ += f"\nPartial function\nFixed args: {all_args}"
        wrapped.__name__ = func.__name__
        wrapped.__module__ = func.__module__

        return wrapped
    return decorator


def vectorize(product = True) -> Callable:
    def decorator(func : Callable) -> Callable:
        signature = inspect.signature(func)
        def wrapped(*args, **kwargs):
            args_to_kwargs = {param : arg for param, arg in zip(signature.parameters, args)}
            all_args = dict(**args_to_kwargs, **kwargs)
            vector_args = {k : v for k, v in all_args.items() if isinstance(v, (list, np.ndarray))}
            shape = [len(v) for v in vector_args.values()]
            scalar_args = {k : v for k, v in all_args.items() if k not in vector_args}

            if product:
                iteration = itertools.product(*vector_args.values())
            else:
                iteration = zip(*vector_args.values())

            results = []
            for it in iteration:
                iteration_args = {k : v for k, v  in zip(vector_args, it)}
                iteration_args.update(scalar_args)
                results.append(func(**iteration_args))
            results = np.reshape(results, shape)
            if not shape:
                results = results.item()
            return results
        wrapped.__signature__ = signature
        wrapped.__annotations__ = generate_annotation(signature)
        wrapped.__doc__ = generate_docstring(signature, return_name=func.__name__)
        wrapped.__name__ = func.__name__
        wrapped.__module__ = func.__module__
        return wrapped
    return decorator


def dummy(func):
    def wrapped(*args, **kwargs):
        return func(*args, **kwargs)
    signature = inspect.signature(func)
    wrapped.__signature__ = signature
    wrapped.__annotations__ = func.__annotations__
    wrapped.__doc__ = func.__doc__
    wrapped.__name__ = func.__name__
    wrapped.__module__ = func.__module__
    return wrapped