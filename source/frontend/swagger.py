from werkzeug.wrappers import Response
import json
from validate import *

def convertType(validator):
    if validator == str :
        return "string"
    elif validator == datetime:
        return "datetime"
    elif validator == time:
        return "time"
    elif validator == int:
        return "integer"
    elif validator == float:
        return "float"
    elif validator == boolean:
        return "boolean"
    else:
        return "string"

def convertBool(boolean):
    if(boolean):
        return "true"
    return "false"

def api_doc(apis, api = None) : 
    response = {}
    response['apiVersion'] = "0.2"
    response['swaggerVersion'] = "1.1"
    response['basePath'] = "http://127.0.0.1:8088"
    response['apis'] = []
    
    if(api) :
        if(api in apis) :
            params = []
            for key, val in apis[api]['arguments'].iteritems():
                param = {}
                param['name'] = key
                param['paramType'] = 'path'
                param['description'] = val.description
                param['dataType'] = convertType(val.validator)
                param['required'] = convertBool(val.required)
                param['allowMultiple'] = convertBool(val.repeated)
                params.append(param)

            response['resourcePath'] = "/"+api
            response['apis'].append({
                    "path" : "/"+api+".{format}",
                    "description" : "",
                    "operations" : [{
                            "httpMethod" : "GET",
                            "summary" : "",
                            "nickname" : api,
                            "responseClass" : "void",
                            "parameters" : params
                            }
                            ]
                                    })

    else:
        for key, val in apis.iteritems() :
            response['apis'].append({"path":"/doc.{format}/"+key, "description" : ""})

    r = Response(json.dumps(response, ensure_ascii=False), mimetype='application/json')
    r.headers.add('Access-Control-Allow-Origin', '*')
    return r

